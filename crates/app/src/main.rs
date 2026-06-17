//! `abclone-verify` — the egui desktop GUI (`docs/03_ARCHITECTURE.md` §6).
//!
//! A thin front-end over `abclone-core`: a project/library sidebar, the four load
//! buttons, Setup / Gate 1 / Gate 2 tabs, colored verdict tables, and a detail
//! pane. All verification logic lives in `abclone-core`; this file is only UI.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use abclone_core::gate1::Gate1Context;
use abclone_core::gate2::{Gate2Context, SeqMode};
use abclone_core::model::{ChainClass, Library, Project};
use abclone_core::seqio::{self, fasta, GroundTruthRow, SeqRecord};
use abclone_core::verdict::{Gate1Verdict, Gate2Verdict};
use abclone_core::{naming, report, workflow};
use eframe::egui;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1100.0, 720.0])
            .with_title("abclone-verify"),
        ..Default::default()
    };
    eframe::run_native(
        "abclone-verify",
        options,
        Box::new(|_cc| {
            // Orders normally carry overhangs; start with that assumption.
            let app = App {
                has_overhangs: true,
                ..Default::default()
            };
            Ok(Box::new(app))
        }),
    )
}

#[derive(PartialEq, Eq, Clone, Copy, Default)]
enum Tab {
    #[default]
    Setup,
    Gate1,
    Gate2,
}

#[derive(Default)]
struct App {
    library: Option<Library>,
    library_path: Option<PathBuf>,
    heavy_vector: Option<String>,

    ground_truth: Vec<GroundTruthRow>,
    order: Vec<SeqRecord>,
    reads: Vec<SeqRecord>,

    order_path: Option<PathBuf>,
    gt_path: Option<PathBuf>,
    reads_path: Option<PathBuf>,

    has_overhangs: bool,
    partial_sanger: bool,

    gate1: Vec<Gate1Verdict>,
    gate2: Vec<Gate2Verdict>,
    sel_g1: Option<usize>,
    sel_g2: Option<usize>,

    tab: Tab,
    status: String,
}

impl App {
    fn project(&self) -> Option<Project> {
        let lib = self.library.as_ref()?;
        Some(workflow::ad_hoc_project(
            lib,
            self.heavy_vector.as_deref(),
            None,
        ))
    }

    fn load_library_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("library", &["json5", "json"])
            .pick_file()
        {
            match Library::load(&path) {
                Ok(lib) => {
                    self.heavy_vector = lib
                        .vectors
                        .iter()
                        .find(|v| v.chain_class == ChainClass::Heavy)
                        .map(|v| v.id.clone());
                    self.status = format!(
                        "Loaded library: {} vector(s), {} overhang set(s)",
                        lib.vectors.len(),
                        lib.overhang_sets.len()
                    );
                    self.library = Some(lib);
                    self.library_path = Some(path);
                }
                Err(e) => self.status = format!("Library load failed: {e}"),
            }
        }
    }

    /// Load the starter library bundled into the binary at build time
    /// (the pre-populated French IgG1 vector + overhang set).
    fn load_bundled_library(&mut self) {
        const BUNDLED: &str = include_str!("../../../reference/example_library.json5");
        match Library::from_json5(BUNDLED) {
            Ok(lib) => {
                self.heavy_vector = lib
                    .vectors
                    .iter()
                    .find(|v| v.chain_class == ChainClass::Heavy)
                    .map(|v| v.id.clone());
                self.status = format!(
                    "Loaded bundled library: {} vector(s), {} overhang set(s)",
                    lib.vectors.len(),
                    lib.overhang_sets.len()
                );
                self.library = Some(lib);
                self.library_path = None;
            }
            Err(e) => self.status = format!("Bundled library failed to parse: {e}"),
        }
    }

    fn load_ground_truth_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("panel", &["csv", "xlsx", "fasta", "fa"])
            .pick_file()
        {
            match load_ground_truth(&path) {
                Ok(gt) => {
                    self.status = format!("Loaded ground truth: {} antibodies", gt.len());
                    self.ground_truth = gt;
                    self.gt_path = Some(path);
                }
                Err(e) => self.status = format!("Ground-truth load failed: {e}"),
            }
        }
    }

    fn load_order_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("order", &["fasta", "fa", "xlsx", "csv"])
            .pick_file()
        {
            match seqio::read_id_seq_auto(&path, "") {
                Ok(recs) => {
                    self.status = format!("Loaded order: {} records", recs.len());
                    self.order = recs;
                    self.order_path = Some(path);
                }
                Err(e) => self.status = format!("Order load failed: {e}"),
            }
        }
    }

    fn load_reads_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("reads", &["fasta", "fa", "ab1"])
            .pick_file()
        {
            let res = if path.extension().and_then(|e| e.to_str()) == Some("ab1") {
                abclone_core::seqio::ab1::read_path(&path)
                    .map(|a| {
                        vec![SeqRecord {
                            id: path
                                .file_stem()
                                .and_then(|s| s.to_str())
                                .unwrap_or("read")
                                .to_string(),
                            sequence: a.bases,
                        }]
                    })
                    .map_err(|e| e.to_string())
            } else {
                fasta::read_path(&path).map_err(|e| e.to_string())
            };
            match res {
                Ok(recs) => {
                    self.status = format!("Loaded sequencing: {} read(s)", recs.len());
                    self.reads = recs;
                    self.reads_path = Some(path);
                }
                Err(e) => self.status = format!("Reads load failed: {e}"),
            }
        }
    }

    fn run_gate1(&mut self) {
        let Some(lib) = self.library.clone() else {
            self.status = "Load a library first".into();
            return;
        };
        let Some(project) = self.project() else {
            return;
        };
        let Some(set) = lib.overhang_set(&project.overhang_set).cloned() else {
            self.status = "Library has no overhang set".into();
            return;
        };
        let ctx = Gate1Context::new(&lib, &project, &set, &self.ground_truth, self.has_overhangs);
        self.gate1 = ctx.run(&self.order);
        let roll = report::rollup_gate1(&self.gate1);
        let pass = roll.iter().filter(|r| r.passed).count();
        self.status = format!("Gate 1: {} of {} antibodies pass", pass, roll.len());
        self.sel_g1 = None;
    }

    fn run_gate2(&mut self) {
        let Some(lib) = self.library.clone() else {
            self.status = "Load a library first".into();
            return;
        };
        let Some(project) = self.project() else {
            return;
        };
        let Some(set) = lib.overhang_set(&project.overhang_set).cloned() else {
            self.status = "Library has no overhang set".into();
            return;
        };
        let cores = workflow::order_cores(&self.order, &lib, &project, self.has_overhangs);
        let mode = if self.partial_sanger {
            SeqMode::PartialSanger
        } else {
            SeqMode::FullPlasmid
        };
        let ctx = Gate2Context::new(&lib, &project, &set, &self.ground_truth, cores, mode);
        self.gate2 = ctx.run(&self.reads);
        let pass = self.gate2.iter().filter(|v| v.passed()).count();
        self.status = format!("Gate 2: {} of {} reads pass", pass, self.gate2.len());
        self.sel_g2 = None;
    }

    fn export_dialog(&mut self) {
        let Some(project) = self.project() else {
            return;
        };
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name("report.html")
            .add_filter("html", &["html"])
            .save_file()
        {
            let html = report::html_report(&project.name, &self.gate1, &self.gate2);
            match std::fs::write(&path, html) {
                Ok(_) => self.status = format!("Wrote {}", path.display()),
                Err(e) => self.status = format!("Export failed: {e}"),
            }
        }
    }
}

impl eframe::App for App {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label("Status:");
                ui.label(&self.status);
            });
        });

        egui::SidePanel::left("sidebar")
            .resizable(true)
            .default_width(280.0)
            .show(ctx, |ui| self.sidebar(ui));

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.tab, Tab::Setup, "Setup");
                ui.selectable_value(&mut self.tab, Tab::Gate1, "Gate 1 · Order QC");
                ui.selectable_value(&mut self.tab, Tab::Gate2, "Gate 2 · Sequencing QC");
            });
            ui.separator();
            match self.tab {
                Tab::Setup => self.setup_tab(ui),
                Tab::Gate1 => self.gate1_tab(ui),
                Tab::Gate2 => self.gate2_tab(ui),
            }
        });
    }
}

impl App {
    fn sidebar(&mut self, ui: &mut egui::Ui) {
        ui.heading("abclone-verify");
        ui.label("Antibody cloning verifier");
        ui.separator();

        ui.label(egui::RichText::new("Library (lab-global)").strong());
        if ui.button("📂 Load Library…").clicked() {
            self.load_library_dialog();
        }
        if ui.button("✨ Use bundled example library").clicked() {
            self.load_bundled_library();
        }
        if let Some(lib) = &self.library {
            ui.label(format!("{} vector(s)", lib.vectors.len()));
            egui::ComboBox::from_label("Heavy vector")
                .selected_text(self.heavy_vector.clone().unwrap_or_else(|| "—".into()))
                .show_ui(ui, |ui| {
                    for v in &lib.vectors {
                        ui.selectable_value(
                            &mut self.heavy_vector,
                            Some(v.id.clone()),
                            &v.display_name,
                        );
                    }
                });
        } else {
            ui.label(egui::RichText::new("none loaded").italics().weak());
        }

        ui.separator();
        ui.label(egui::RichText::new("Project inputs").strong());
        input_row(
            ui,
            "Ground truth",
            self.gt_path.as_deref(),
            self.ground_truth.len(),
        );
        input_row(
            ui,
            "IDT order",
            self.order_path.as_deref(),
            self.order.len(),
        );
        input_row(
            ui,
            "Sequencing",
            self.reads_path.as_deref(),
            self.reads.len(),
        );

        ui.separator();
        ui.checkbox(&mut self.has_overhangs, "Order has overhangs");
        ui.checkbox(&mut self.partial_sanger, "Partial Sanger reads");

        ui.separator();
        if ui.button("💾 Export HTML report…").clicked() {
            self.export_dialog();
        }
    }

    fn setup_tab(&mut self, ui: &mut egui::Ui) {
        ui.label("Load the four inputs, then run a gate. The library is lab-global; the others are per-project.");
        ui.add_space(6.0);
        ui.horizontal_wrapped(|ui| {
            if ui.button("📂 Load Library").clicked() {
                self.load_library_dialog();
            }
            if ui.button("📂 Load Ground Truth").clicked() {
                self.load_ground_truth_dialog();
            }
            if ui.button("📂 Load IDT Order").clicked() {
                self.load_order_dialog();
            }
            if ui.button("📂 Load Sequencing Results").clicked() {
                self.load_reads_dialog();
            }
        });
        ui.separator();

        ui.label(egui::RichText::new("Parsed order records").strong());
        let profile = self
            .library
            .as_ref()
            .and_then(|l| l.naming_profiles.first().cloned())
            .unwrap_or_else(naming::default_profile);
        egui::ScrollArea::vertical()
            .max_height(420.0)
            .show(ui, |ui| {
                egui::Grid::new("records")
                    .striped(true)
                    .num_columns(4)
                    .show(ui, |ui| {
                        ui.label(egui::RichText::new("Record").strong());
                        ui.label(egui::RichText::new("Antibody").strong());
                        ui.label(egui::RichText::new("Chain").strong());
                        ui.label(egui::RichText::new("Confidence").strong());
                        ui.end_row();
                        for r in &self.order {
                            let p = naming::parse_name(&r.id, &profile);
                            ui.label(&r.id);
                            ui.label(&p.ab_id);
                            ui.label(p.chain_class.as_str());
                            let c = if p.needs_confirmation {
                                "⚠ confirm"
                            } else {
                                "ok"
                            };
                            ui.label(c);
                            ui.end_row();
                        }
                    });
            });
    }

    fn gate1_tab(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if ui.button("▶ Run Gate 1").clicked() {
                self.run_gate1();
            }
            ui.label(format!("{} verdicts", self.gate1.len()));
        });
        ui.separator();
        let mut clicked: Option<usize> = None;
        egui::ScrollArea::vertical()
            .max_height(360.0)
            .show(ui, |ui| {
                egui::Grid::new("g1")
                    .striped(true)
                    .num_columns(5)
                    .show(ui, |ui| {
                        for s in ["Record", "Antibody", "Chain", "Verdict", "Reads through"] {
                            ui.label(egui::RichText::new(s).strong());
                        }
                        ui.end_row();
                        for (i, v) in self.gate1.iter().enumerate() {
                            if ui
                                .selectable_label(self.sel_g1 == Some(i), &v.record_id)
                                .clicked()
                            {
                                clicked = Some(i);
                            }
                            ui.label(&v.ab_id);
                            ui.label(&v.chain_class);
                            ui.colored_label(verdict_color(v.kind.label()), v.kind.label());
                            ui.label(
                                v.reads_through
                                    .map(|b| b.to_string())
                                    .unwrap_or_else(|| "—".into()),
                            );
                            ui.end_row();
                        }
                    });
            });
        if clicked.is_some() {
            self.sel_g1 = clicked;
        }
        if let Some(i) = self.sel_g1 {
            if let Some(v) = self.gate1.get(i) {
                ui.separator();
                ui.label(egui::RichText::new("Detail").strong());
                ui.label(format!("{} — {}", v.record_id, v.kind.label()));
                ui.label(&v.reason);
                if let Some(c) = v.core_len {
                    ui.label(format!("core length: {c} nt"));
                }
                if let Some(s) = v.premature_stop_aa {
                    ui.label(format!("premature stop at aa {s}"));
                }
            }
        }
    }

    fn gate2_tab(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if ui.button("▶ Run Gate 2").clicked() {
                self.run_gate2();
            }
            ui.label(format!("{} verdicts", self.gate2.len()));
        });
        ui.separator();
        let mut clicked: Option<usize> = None;
        egui::ScrollArea::vertical()
            .max_height(360.0)
            .show(ui, |ui| {
                egui::Grid::new("g2")
                    .striped(true)
                    .num_columns(5)
                    .show(ui, |ui| {
                        for s in ["Record", "Verdict", "Backbone", "Identity", "Reads through"] {
                            ui.label(egui::RichText::new(s).strong());
                        }
                        ui.end_row();
                        for (i, v) in self.gate2.iter().enumerate() {
                            if ui
                                .selectable_label(self.sel_g2 == Some(i), &v.record_id)
                                .clicked()
                            {
                                clicked = Some(i);
                            }
                            ui.colored_label(verdict_color(v.kind.label()), v.kind.label());
                            ui.label(v.backbone_vector.clone().unwrap_or_else(|| "—".into()));
                            ui.label(
                                v.backbone_identity
                                    .map(|x| format!("{:.1}%", x * 100.0))
                                    .unwrap_or_else(|| "—".into()),
                            );
                            ui.label(
                                v.reads_through
                                    .map(|b| b.to_string())
                                    .unwrap_or_else(|| "—".into()),
                            );
                            ui.end_row();
                        }
                    });
            });
        if clicked.is_some() {
            self.sel_g2 = clicked;
        }
        if let Some(i) = self.sel_g2 {
            if let Some(v) = self.gate2.get(i) {
                ui.separator();
                ui.label(egui::RichText::new("Detail").strong());
                ui.label(format!("{} — {}", v.record_id, v.kind.label()));
                ui.label(&v.reason);
                if !v.mutations.is_empty() {
                    let muts: Vec<String> = v
                        .mutations
                        .iter()
                        .map(|m| format!("{}{}→{}", m.wt, m.position_aa, m.mut_aa))
                        .collect();
                    ui.label(format!("mutations: {}", muts.join(", ")));
                }
                if let Some(s) = &v.suspected_identity {
                    ui.label(format!("suspected correct identity: {s}"));
                }
            }
        }
    }
}

fn input_row(ui: &mut egui::Ui, label: &str, path: Option<&Path>, n: usize) {
    ui.horizontal(|ui| {
        let mark = if path.is_some() { "✔" } else { "—" };
        ui.label(format!("{mark} {label}"));
        if path.is_some() {
            ui.weak(format!("({n})"));
        }
    });
}

fn verdict_color(label: &str) -> egui::Color32 {
    match label {
        "PASS" => egui::Color32::from_rgb(26, 127, 55),
        "NO_GROUND_TRUTH" | "SILENT_VARIANT" | "GC_WARNING" | "RARE_CODON_WARNING" => {
            egui::Color32::from_rgb(154, 103, 0)
        }
        _ => egui::Color32::from_rgb(207, 34, 46),
    }
}

fn load_ground_truth(path: &Path) -> Result<Vec<GroundTruthRow>, String> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase());
    if matches!(ext.as_deref(), Some("fasta") | Some("fa")) {
        let recs = fasta::read_path(path).map_err(|e| e.to_string())?;
        let prof = naming::default_profile();
        let mut map: BTreeMap<String, GroundTruthRow> = BTreeMap::new();
        for r in recs {
            let np = naming::parse_name(&r.id, &prof);
            let row = map
                .entry(np.ab_id.clone())
                .or_insert_with(|| GroundTruthRow {
                    ab_id: np.ab_id.clone(),
                    heavy: None,
                    light: None,
                });
            match np.chain_class {
                ChainClass::Heavy => row.heavy = Some(r.sequence),
                _ => row.light = Some(r.sequence),
            }
        }
        Ok(map.into_values().collect())
    } else {
        seqio::tabular::read_ground_truth_table(path, &BTreeMap::new()).map_err(|e| e.to_string())
    }
}
