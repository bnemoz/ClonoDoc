//! `clonodoc` — the egui desktop GUI (`docs/03_ARCHITECTURE.md` §6).
//!
//! A thin front-end over `clonodoc-core`: a library/project sidebar, the load
//! buttons, a guided library builder, Setup / Gate 1 / Gate 2 tabs, colored
//! verdict tables, and a detail pane. All verification logic lives in
//! `clonodoc-core`; this file is only UI.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use clonodoc_core::gate1::Gate1Context;
use clonodoc_core::gate2::{Gate2Context, SeqMode};
use clonodoc_core::model::{ChainClass, Library, Locus, Project};
use clonodoc_core::seqio::genbank::GbRecord;
use clonodoc_core::seqio::{self, fasta, GroundTruthRow, SeqRecord};
use clonodoc_core::verdict::{Gate1Verdict, Gate2Verdict};
use clonodoc_core::{naming, report, workflow};
use eframe::egui;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1180.0, 760.0])
            .with_title("ClonoDoc"),
        ..Default::default()
    };
    eframe::run_native(
        "clonodoc",
        options,
        Box::new(|_cc| Ok(Box::new(App::new()))),
    )
}

#[derive(PartialEq, Eq, Clone, Copy, Default)]
enum Tab {
    #[default]
    Setup,
    Gate1,
    Gate2,
}

/// In-progress "add vector from GenBank" form state.
#[derive(Default)]
struct AddVectorForm {
    open: bool,
    gb_path: Option<PathBuf>,
    gb: Option<GbRecord>,
    id: String,
    display: String,
    isotype: String,
    locus: Locus,
    overhang_set: String,
    error: String,
}

#[derive(Default)]
struct App {
    library: Option<Library>,
    library_path: Option<PathBuf>,
    library_dirty: bool,
    heavy_vector: Option<String>,

    ground_truth: Vec<GroundTruthRow>,
    order: Vec<SeqRecord>,
    reads: Vec<SeqRecord>,

    order_path: Option<PathBuf>,
    // Ground-truth source + column mapping (for custom-named columns).
    gt_path: Option<PathBuf>,
    gt_is_fasta: bool,
    gt_headers: Vec<String>,
    gt_ab_col: Option<String>,
    gt_heavy_col: Option<String>,
    gt_light_col: Option<String>,
    // Sequencing source: a set of explicit files and/or a folder (re-read on Run).
    reads_files: Vec<PathBuf>,
    reads_dir: Option<PathBuf>,

    has_overhangs: bool,
    partial_sanger: bool,

    add_vector: AddVectorForm,

    gate1: Vec<Gate1Verdict>,
    gate2: Vec<Gate2Verdict>,
    sel_g1: Option<usize>,
    sel_g2: Option<usize>,

    tab: Tab,
    status: String,
}

impl App {
    fn new() -> Self {
        App {
            has_overhangs: true,
            status: "Load or build a library to begin. New here? Library ▸ New library, then Add vector from GenBank.".into(),
            ..Default::default()
        }
    }

    fn project(&self) -> Option<Project> {
        let lib = self.library.as_ref()?;
        Some(workflow::ad_hoc_project(
            lib,
            self.heavy_vector.as_deref(),
            None,
        ))
    }

    // ---- Library ---------------------------------------------------------

    fn set_library(&mut self, lib: Library, path: Option<PathBuf>) {
        self.heavy_vector = lib
            .vectors
            .iter()
            .find(|v| v.chain_class == ChainClass::Heavy)
            .map(|v| v.id.clone());
        self.status = format!(
            "Library: {} vector(s), {} overhang set(s)",
            lib.vectors.len(),
            lib.overhang_sets.len()
        );
        self.library_dirty = false;
        self.library = Some(lib);
        self.library_path = path;
    }

    fn load_library_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("library", &["json5", "json"])
            .pick_file()
        {
            match Library::load(&path) {
                Ok(lib) => self.set_library(lib, Some(path)),
                Err(e) => self.status = format!("Library load failed: {e}"),
            }
        }
    }

    /// The starter library bundled into the binary (French IgG1 vector + overhangs).
    fn load_bundled_library(&mut self) {
        const BUNDLED: &str = include_str!("../../../reference/example_library.json5");
        match Library::from_json5(BUNDLED) {
            Ok(lib) => self.set_library(lib, None),
            Err(e) => self.status = format!("Bundled library failed to parse: {e}"),
        }
    }

    /// A new library seeded with the bundled naming profile + overhang set(s) but
    /// no vectors — the starting point for building your own from scratch.
    fn new_library(&mut self) {
        const BUNDLED: &str = include_str!("../../../reference/example_library.json5");
        let mut lib = Library::from_json5(BUNDLED).unwrap_or_else(|_| Library::empty());
        lib.vectors.clear();
        self.set_library(lib, None);
        self.library_dirty = true;
        self.status =
            "New library created (overhang set + naming profile kept; no vectors). Add a vector from GenBank.".into();
    }

    fn save_library_dialog(&mut self) {
        let Some(lib) = &self.library else {
            self.status = "No library to save".into();
            return;
        };
        let mut dialog = rfd::FileDialog::new().add_filter("library", &["json5", "json"]);
        if let Some(p) = &self.library_path {
            dialog = dialog.set_file_name(
                p.file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or("library.json5"),
            );
        } else {
            dialog = dialog.set_file_name("library.json5");
        }
        if let Some(path) = dialog.save_file() {
            match lib.save(&path) {
                Ok(_) => {
                    self.status = format!("Saved library to {}", path.display());
                    self.library_path = Some(path);
                    self.library_dirty = false;
                }
                Err(e) => self.status = format!("Save failed: {e}"),
            }
        }
    }

    fn open_add_vector(&mut self) {
        let default_set = self
            .library
            .as_ref()
            .and_then(|l| l.overhang_sets.first())
            .map(|o| o.id.clone())
            .unwrap_or_default();
        self.add_vector = AddVectorForm {
            open: true,
            isotype: "IgG1".into(),
            overhang_set: default_set,
            ..Default::default()
        };
    }

    fn add_vector_pick_gb(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("GenBank", &["gb", "gbk", "genbank"])
            .pick_file()
        {
            match seqio::genbank::read_path(&path) {
                Ok(rec) => {
                    let stem = path
                        .file_stem()
                        .and_then(|s| s.to_str())
                        .unwrap_or("vector")
                        .to_string();
                    if self.add_vector.id.is_empty() {
                        self.add_vector.id = sanitize_id(&stem);
                    }
                    if self.add_vector.display.is_empty() {
                        self.add_vector.display = stem;
                    }
                    self.add_vector.gb = Some(rec);
                    self.add_vector.gb_path = Some(path);
                    self.add_vector.error.clear();
                }
                Err(e) => self.add_vector.error = format!("Could not parse GenBank: {e}"),
            }
        }
    }

    fn commit_add_vector(&mut self) {
        let form = &self.add_vector;
        let Some(lib) = self.library.as_mut() else {
            self.add_vector.error = "Create or load a library first".into();
            return;
        };
        let Some(gb) = &form.gb else {
            self.add_vector.error = "Pick a GenBank (.gb) file first".into();
            return;
        };
        if form.id.trim().is_empty() {
            self.add_vector.error = "Give the vector an id".into();
            return;
        }
        if lib.vectors.iter().any(|v| v.id == form.id) {
            self.add_vector.error = format!("A vector with id '{}' already exists", form.id);
            return;
        }
        let Some(set) = lib.overhang_set(&form.overhang_set).cloned() else {
            self.add_vector.error = "Select an overhang set (the library has none)".into();
            return;
        };
        let v = workflow::vector_from_genbank(
            gb,
            form.id.trim(),
            if form.display.trim().is_empty() {
                form.id.trim()
            } else {
                form.display.trim()
            },
            form.locus.chain_class(),
            form.isotype.trim(),
            &set,
            form.locus,
            "gui",
        );
        let summary = format!(
            "Added vector '{}' ({}, {} bp): insertion site {}/{}, anchor {}",
            v.id,
            v.chain_class.as_str(),
            v.length,
            v.insertion_site.oh5_end,
            v.insertion_site.oh3_start,
            v.constant_anchor_aa
        );
        if v.chain_class == ChainClass::Heavy {
            self.heavy_vector = Some(v.id.clone());
        }
        lib.vectors.push(v);
        self.library_dirty = true;
        self.add_vector.open = false;
        self.status = format!("{summary} — remember to Save library.");
    }

    // ---- Inputs ----------------------------------------------------------

    fn load_order_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("order", &["fasta", "fa", "fna", "xlsx", "csv"])
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

    fn load_ground_truth_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new()
            .add_filter("panel", &["csv", "xlsx", "fasta", "fa"])
            .pick_file()
        {
            let is_fasta = matches!(
                path.extension()
                    .and_then(|e| e.to_str())
                    .map(|e| e.to_ascii_lowercase())
                    .as_deref(),
                Some("fasta") | Some("fa") | Some("fna")
            );
            self.gt_path = Some(path.clone());
            self.gt_is_fasta = is_fasta;
            if is_fasta {
                self.gt_headers.clear();
                self.parse_ground_truth();
            } else {
                match seqio::tabular::read_headers(&path) {
                    Ok(headers) => {
                        // Auto-guess the mapping from the headers, then parse.
                        self.gt_ab_col = guess_col(
                            &headers,
                            &[
                                "ab_id",
                                "abid",
                                "antibody",
                                "id",
                                "name",
                                "clone",
                                "sequence_id",
                            ],
                        );
                        self.gt_heavy_col = guess_col(
                            &headers,
                            &["h_seq", "heavy", "heavychain", "hc", "vh", "heavy_seq"],
                        );
                        self.gt_light_col = guess_col(
                            &headers,
                            &["l_seq", "light", "lightchain", "lc", "vl", "light_seq"],
                        );
                        self.gt_headers = headers;
                        self.parse_ground_truth();
                    }
                    Err(e) => self.status = format!("Could not read columns: {e}"),
                }
            }
        }
    }

    /// (Re)parse the ground truth from its stored path + current column mapping.
    fn parse_ground_truth(&mut self) {
        let Some(path) = self.gt_path.clone() else {
            return;
        };
        if self.gt_is_fasta {
            match load_ground_truth_fasta(&path) {
                Ok(gt) => {
                    self.ground_truth = gt;
                    self.status = format!(
                        "Ground truth: {} antibodies (FASTA)",
                        self.ground_truth.len()
                    );
                }
                Err(e) => self.status = format!("Ground-truth load failed: {e}"),
            }
            return;
        }
        let mut overrides: BTreeMap<String, String> = BTreeMap::new();
        if let Some(c) = &self.gt_ab_col {
            overrides.insert("ab_id".into(), c.clone());
        }
        if let Some(c) = &self.gt_heavy_col {
            overrides.insert("heavy".into(), c.clone());
        }
        if let Some(c) = &self.gt_light_col {
            overrides.insert("light".into(), c.clone());
        }
        match seqio::tabular::read_ground_truth_table(&path, &overrides) {
            Ok(gt) => {
                let with_seq = gt
                    .iter()
                    .filter(|r| r.heavy.is_some() || r.light.is_some())
                    .count();
                self.ground_truth = gt;
                self.status = format!(
                    "Ground truth: {} antibodies ({} with a heavy/light sequence)",
                    self.ground_truth.len(),
                    with_seq
                );
                if with_seq == 0 {
                    self.status
                        .push_str(" — check the column mapping (no sequences matched).");
                }
            }
            Err(e) => self.status = format!("Ground-truth parse failed: {e}"),
        }
    }

    fn load_reads_files_dialog(&mut self) {
        if let Some(paths) = rfd::FileDialog::new()
            .add_filter("reads", &["fasta", "fa", "fna", "ab1"])
            .pick_files()
        {
            self.reads_files = paths;
            self.reads_dir = None;
            self.reload_reads();
        }
    }

    fn load_reads_folder_dialog(&mut self) {
        if let Some(dir) = rfd::FileDialog::new().pick_folder() {
            self.reads_dir = Some(dir);
            self.reads_files.clear();
            self.reload_reads();
        }
    }

    /// Re-read all sequencing reads from the current files/folder.
    fn reload_reads(&mut self) {
        let mut all = Vec::new();
        let mut errors = Vec::new();
        let paths: Vec<PathBuf> = if let Some(dir) = &self.reads_dir {
            list_read_files(dir)
        } else {
            self.reads_files.clone()
        };
        let n_files = paths.len();
        for p in &paths {
            match read_one_reads_file(p) {
                Ok(mut recs) => all.append(&mut recs),
                Err(e) => errors.push(format!("{}: {e}", p.display())),
            }
        }
        self.reads = all;
        self.status = if errors.is_empty() {
            format!(
                "Loaded sequencing: {} read(s) from {} file(s)",
                self.reads.len(),
                n_files
            )
        } else {
            format!(
                "Loaded {} read(s) from {} file(s); {} file(s) failed",
                self.reads.len(),
                n_files,
                errors.len()
            )
        };
    }

    /// Re-read order + ground truth + reads from disk so that fixing a file and
    /// re-running a gate always reflects the latest content.
    fn reload_inputs(&mut self) {
        if let Some(p) = self.order_path.clone() {
            if let Ok(recs) = seqio::read_id_seq_auto(&p, "") {
                self.order = recs;
            }
        }
        if self.gt_path.is_some() {
            self.parse_ground_truth();
        }
        if self.reads_dir.is_some() || !self.reads_files.is_empty() {
            self.reload_reads();
        }
    }

    // ---- Gates -----------------------------------------------------------

    fn run_gate1(&mut self) {
        self.reload_inputs();
        let Some(lib) = self.library.clone() else {
            self.status = "Load or build a library first".into();
            return;
        };
        let Some(project) = self.project() else {
            return;
        };
        let Some(set) = lib.overhang_set(&project.overhang_set).cloned() else {
            self.status = "Library has no overhang set".into();
            return;
        };
        if self.order.is_empty() {
            self.status = "Load an IDT order first".into();
            return;
        }
        let ctx = Gate1Context::new(&lib, &project, &set, &self.ground_truth, self.has_overhangs);
        self.gate1 = ctx.run(&self.order);
        let roll = report::rollup_gate1(&self.gate1);
        let pass = roll.iter().filter(|r| r.passed).count();
        self.status = format!(
            "Gate 1: {} of {} antibodies pass all chains",
            pass,
            roll.len()
        );
        self.sel_g1 = None;
        self.tab = Tab::Gate1;
    }

    fn run_gate2(&mut self) {
        self.reload_inputs();
        let Some(lib) = self.library.clone() else {
            self.status = "Load or build a library first".into();
            return;
        };
        let Some(project) = self.project() else {
            return;
        };
        let Some(set) = lib.overhang_set(&project.overhang_set).cloned() else {
            self.status = "Library has no overhang set".into();
            return;
        };
        if self.reads.is_empty() {
            self.status = "Load sequencing reads first".into();
            return;
        }
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
        self.tab = Tab::Gate2;
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
            ui.horizontal_wrapped(|ui| {
                ui.label(egui::RichText::new("Status:").strong());
                ui.label(&self.status);
            });
        });

        egui::SidePanel::left("sidebar")
            .resizable(true)
            .default_width(300.0)
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

        self.add_vector_window(ctx);
    }
}

impl App {
    fn sidebar(&mut self, ui: &mut egui::Ui) {
        ui.heading("ClonoDoc");
        ui.label("Antibody cloning verifier");
        ui.separator();

        ui.label(egui::RichText::new("Library (lab-global)").strong());
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("📂 Load…")
                .on_hover_text("Open an existing library.json5")
                .clicked()
            {
                self.load_library_dialog();
            }
            if ui
                .button("🆕 New")
                .on_hover_text("Start a new library (keeps overhangs + naming, no vectors)")
                .clicked()
            {
                self.new_library();
            }
            if ui
                .button("✨ Example")
                .on_hover_text("Load the bundled French IgG1 example library")
                .clicked()
            {
                self.load_bundled_library();
            }
        });
        let has_lib = self.library.is_some();
        ui.add_enabled_ui(has_lib, |ui| {
            ui.horizontal_wrapped(|ui| {
                if ui.button("➕ Add vector from GenBank…").clicked() {
                    self.open_add_vector();
                }
                let save_label = if self.library_dirty {
                    "💾 Save library* …"
                } else {
                    "💾 Save library…"
                };
                if ui.button(save_label).clicked() {
                    self.save_library_dialog();
                }
            });
        });

        if let Some(lib) = &self.library {
            ui.add_space(4.0);
            ui.label(format!(
                "{} vector(s), {} overhang set(s)",
                lib.vectors.len(),
                lib.overhang_sets.len()
            ));
            if lib.vectors.is_empty() {
                ui.label(
                    egui::RichText::new("No vectors yet — use “Add vector from GenBank…”.")
                        .italics()
                        .weak(),
                );
            } else {
                let combo = egui::ComboBox::from_id_salt("heavy_vec").selected_text(
                    self.heavy_vector
                        .as_deref()
                        .and_then(|id| lib.vector(id))
                        .map(|v| v.display_name.clone())
                        .unwrap_or_else(|| "— none —".into()),
                );
                ui.horizontal(|ui| {
                    ui.label("Heavy-chain backbone:")
                        .on_hover_text("Which library vector is the backbone for HEAVY chains.\nκ/λ light vectors are auto-detected from the read, so they are not selected here.");
                    combo.show_ui(ui, |ui| {
                        for v in lib.vectors.iter().filter(|v| v.chain_class == ChainClass::Heavy) {
                            ui.selectable_value(&mut self.heavy_vector, Some(v.id.clone()), &v.display_name);
                        }
                    });
                });
            }
        } else {
            ui.label(egui::RichText::new("No library loaded").italics().weak());
        }

        ui.separator();
        ui.label(egui::RichText::new("Project inputs").strong());
        input_row(
            ui,
            "Ground truth",
            self.gt_path.is_some(),
            self.ground_truth.len(),
        );
        input_row(ui, "IDT order", self.order_path.is_some(), self.order.len());
        let reads_loaded = self.reads_dir.is_some() || !self.reads_files.is_empty();
        input_row(ui, "Sequencing", reads_loaded, self.reads.len());

        ui.separator();
        ui.checkbox(&mut self.has_overhangs, "Order has overhangs")
            .on_hover_text("Tick if the ordered sequences include the 5′/3′ overhangs (e.g. IDT FASTA). Untick for bare optimized cores.");
        ui.checkbox(&mut self.partial_sanger, "Partial Sanger reads")
            .on_hover_text(
                "Tick for junction/insert-only reads (backbone reported as flanks_only).",
            );

        ui.separator();
        if ui.button("💾 Export HTML report…").clicked() {
            self.export_dialog();
        }
    }

    fn add_vector_window(&mut self, ctx: &egui::Context) {
        if !self.add_vector.open {
            return;
        }
        let mut open = self.add_vector.open;
        let mut do_pick = false;
        let mut do_commit = false;
        let mut do_cancel = false;
        let overhang_ids: Vec<String> = self
            .library
            .as_ref()
            .map(|l| l.overhang_sets.iter().map(|o| o.id.clone()).collect())
            .unwrap_or_default();

        egui::Window::new("Add vector from GenBank")
            .open(&mut open)
            .resizable(true)
            .default_width(460.0)
            .show(ctx, |ui| {
                let f = &mut self.add_vector;
                ui.horizontal(|ui| {
                    if ui.button("📂 Choose .gb file…").clicked() {
                        do_pick = true;
                    }
                    match &f.gb_path {
                        Some(p) => ui.label(p.file_name().and_then(|n| n.to_str()).unwrap_or("")),
                        None => ui.label(egui::RichText::new("no file chosen").weak()),
                    };
                });
                if let Some(gb) = &f.gb {
                    ui.label(format!(
                        "{} — {} bp, {}",
                        gb.name,
                        gb.length,
                        if gb.circular { "circular" } else { "linear" }
                    ));
                }
                ui.separator();
                egui::Grid::new("av_form").num_columns(2).show(ui, |ui| {
                    ui.label("Vector id");
                    ui.text_edit_singleline(&mut f.id);
                    ui.end_row();
                    ui.label("Display name");
                    ui.text_edit_singleline(&mut f.display);
                    ui.end_row();
                    ui.label("Isotype");
                    ui.text_edit_singleline(&mut f.isotype);
                    ui.end_row();
                    ui.label("Locus / chain");
                    egui::ComboBox::from_id_salt("av_locus")
                        .selected_text(locus_label(f.locus))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut f.locus, Locus::Igh, locus_label(Locus::Igh));
                            ui.selectable_value(&mut f.locus, Locus::Igk, locus_label(Locus::Igk));
                            ui.selectable_value(&mut f.locus, Locus::Igl, locus_label(Locus::Igl));
                        });
                    ui.end_row();
                    ui.label("Overhang set");
                    egui::ComboBox::from_id_salt("av_set")
                        .selected_text(if f.overhang_set.is_empty() { "—".into() } else { f.overhang_set.clone() })
                        .show_ui(ui, |ui| {
                            for id in &overhang_ids {
                                ui.selectable_value(&mut f.overhang_set, id.clone(), id);
                            }
                        });
                    ui.end_row();
                });
                ui.label(
                    egui::RichText::new(
                        "The insertion site and constant-region anchor are computed automatically by locating the overhangs in the vector (falling back to the GenBank feature boundaries).",
                    )
                    .weak()
                    .small(),
                );
                if !f.error.is_empty() {
                    ui.colored_label(egui::Color32::from_rgb(207, 34, 46), &f.error);
                }
                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Add vector").clicked() {
                        do_commit = true;
                    }
                    if ui.button("Cancel").clicked() {
                        do_cancel = true;
                    }
                });
            });

        if do_pick {
            self.add_vector_pick_gb();
        }
        if do_commit {
            self.commit_add_vector();
        }
        if do_cancel {
            self.add_vector.open = false;
            return;
        }
        // Honor the window's close [x] unless we already committed/closed.
        if self.add_vector.open {
            self.add_vector.open = open;
        }
    }

    fn setup_tab(&mut self, ui: &mut egui::Ui) {
        ui.label("Load the inputs, then run a gate. The library is lab-global; ground truth / order / sequencing are per-project.");
        ui.add_space(6.0);
        ui.horizontal_wrapped(|ui| {
            if ui.button("📂 Ground Truth…").clicked() {
                self.load_ground_truth_dialog();
            }
            if ui.button("📂 IDT Order…").clicked() {
                self.load_order_dialog();
            }
            if ui.button("📂 Sequencing files…").clicked() {
                self.load_reads_files_dialog();
            }
            if ui.button("📁 Sequencing folder…").clicked() {
                self.load_reads_folder_dialog();
            }
        });

        // Ground-truth column mapping for tables with custom column names.
        if !self.gt_headers.is_empty() {
            ui.separator();
            ui.label(egui::RichText::new("Ground-truth columns").strong());
            ui.label(
                egui::RichText::new(
                    "Map your columns to antibody id + heavy/light sequence, then Apply.",
                )
                .weak()
                .small(),
            );
            let headers = self.gt_headers.clone();
            let mut changed = false;
            changed |= column_picker(ui, "Antibody id", "gt_ab", &headers, &mut self.gt_ab_col);
            changed |= column_picker(
                ui,
                "Heavy sequence",
                "gt_h",
                &headers,
                &mut self.gt_heavy_col,
            );
            changed |= column_picker(
                ui,
                "Light sequence",
                "gt_l",
                &headers,
                &mut self.gt_light_col,
            );
            if ui.button("Apply mapping").clicked() || changed {
                self.parse_ground_truth();
            }
        }

        ui.separator();
        ui.label(egui::RichText::new("Parsed order records").strong());
        let profile = self
            .library
            .as_ref()
            .and_then(|l| l.naming_profiles.first().cloned())
            .unwrap_or_else(naming::default_profile);
        egui::ScrollArea::vertical()
            .max_height(380.0)
            .show(ui, |ui| {
                egui::Grid::new("records")
                    .striped(true)
                    .num_columns(4)
                    .show(ui, |ui| {
                        for s in ["Record", "Antibody", "Chain", "Confidence"] {
                            ui.label(egui::RichText::new(s).strong());
                        }
                        ui.end_row();
                        for r in &self.order {
                            let p = naming::parse_name(&r.id, &profile);
                            ui.label(&r.id);
                            ui.label(&p.ab_id);
                            ui.label(p.chain_class.as_str());
                            ui.label(if p.needs_confirmation {
                                "⚠ confirm"
                            } else {
                                "ok"
                            });
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
            ui.label("(re-reads files from disk on each run)");
            ui.label(format!("· {} verdicts", self.gate1.len()));
        });
        ui.separator();
        let mut clicked = None;
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
                            ui.label(opt_bool(v.reads_through));
                            ui.end_row();
                        }
                    });
            });
        if clicked.is_some() {
            self.sel_g1 = clicked;
        }
        if let Some(v) = self.sel_g1.and_then(|i| self.gate1.get(i)) {
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

    fn gate2_tab(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if ui.button("▶ Run Gate 2").clicked() {
                self.run_gate2();
            }
            ui.label("(re-reads files from disk on each run)");
            ui.label(format!("· {} verdicts", self.gate2.len()));
        });
        ui.separator();
        let mut clicked = None;
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
                            ui.label(opt_bool(v.reads_through));
                            ui.end_row();
                        }
                    });
            });
        if clicked.is_some() {
            self.sel_g2 = clicked;
        }
        if let Some(v) = self.sel_g2.and_then(|i| self.gate2.get(i)) {
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

// --- free helpers -----------------------------------------------------------

fn input_row(ui: &mut egui::Ui, label: &str, loaded: bool, n: usize) {
    ui.horizontal(|ui| {
        ui.label(format!("{} {label}", if loaded { "✔" } else { "—" }));
        if loaded {
            ui.weak(format!("({n})"));
        }
    });
}

fn opt_bool(b: Option<bool>) -> String {
    b.map(|x| x.to_string()).unwrap_or_else(|| "—".into())
}

fn locus_label(l: Locus) -> &'static str {
    match l {
        Locus::Igh => "IGH (heavy)",
        Locus::Igk => "IGK (kappa)",
        Locus::Igl => "IGL (lambda)",
    }
}

fn sanitize_id(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() {
                c.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect()
}

/// A dropdown to map a logical field to one of the file's columns (or none).
fn column_picker(
    ui: &mut egui::Ui,
    label: &str,
    id: &str,
    headers: &[String],
    selected: &mut Option<String>,
) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label(format!("{label}:"));
        egui::ComboBox::from_id_salt(id)
            .selected_text(selected.clone().unwrap_or_else(|| "— none —".into()))
            .show_ui(ui, |ui| {
                if ui
                    .selectable_label(selected.is_none(), "— none —")
                    .clicked()
                {
                    *selected = None;
                    changed = true;
                }
                for h in headers {
                    if ui
                        .selectable_label(selected.as_deref() == Some(h.as_str()), h)
                        .clicked()
                    {
                        *selected = Some(h.clone());
                        changed = true;
                    }
                }
            });
    });
    changed
}

fn guess_col(headers: &[String], aliases: &[&str]) -> Option<String> {
    let norm = |s: &str| -> String {
        s.trim()
            .to_ascii_lowercase()
            .chars()
            .filter(|c| c.is_ascii_alphanumeric())
            .collect()
    };
    for a in aliases {
        let an = norm(a);
        if let Some(h) = headers.iter().find(|h| norm(h) == an) {
            return Some(h.clone());
        }
    }
    None
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

/// Read one sequencing file (AB1 or FASTA) into records.
fn read_one_reads_file(path: &Path) -> Result<Vec<SeqRecord>, String> {
    if path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("ab1"))
        .unwrap_or(false)
    {
        let a = seqio::ab1::read_path(path).map_err(|e| e.to_string())?;
        let id = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("read")
            .to_string();
        Ok(vec![SeqRecord {
            id,
            sequence: a.bases,
        }])
    } else {
        fasta::read_path(path).map_err(|e| e.to_string())
    }
}

/// List the AB1/FASTA files directly inside a folder.
fn list_read_files(dir: &Path) -> Vec<PathBuf> {
    let mut out = Vec::new();
    if let Ok(entries) = std::fs::read_dir(dir) {
        for e in entries.flatten() {
            let p = e.path();
            if let Some(ext) = p.extension().and_then(|x| x.to_str()) {
                if matches!(
                    ext.to_ascii_lowercase().as_str(),
                    "ab1" | "fasta" | "fa" | "fna"
                ) {
                    out.push(p);
                }
            }
        }
    }
    out.sort();
    out
}

fn load_ground_truth_fasta(path: &Path) -> Result<Vec<GroundTruthRow>, String> {
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
}
