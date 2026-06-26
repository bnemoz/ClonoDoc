//! `clonodoc` — the egui desktop GUI, styled to the "Bench" dark design
//! (`design/README.md`). A thin front-end over `clonodoc-core`: a library/project
//! sidebar, four top tabs (Reference Library · Project Files · In-silico Check ·
//! Wetlab Verify), and batch verdict tables that drill into a rich per-construct
//! detail (plasmid map, junction viewer, alignment). All verification logic lives
//! in `clonodoc-core`; this file is only UI + data-driven drawing.
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

mod theme;
mod viz;

use clonodoc_core::assemble;
use clonodoc_core::gate1::Gate1Context;
use clonodoc_core::gate2::{Gate2Context, SeqMode};
use clonodoc_core::model::{ChainClass, Library, Locus, Project, Vector};
use clonodoc_core::seqio::genbank::GbRecord;
use clonodoc_core::seqio::{self, fasta, GroundTruthRow, SeqRecord};
use clonodoc_core::verdict::{Gate1Verdict, Gate2Verdict};
use clonodoc_core::{naming, report, seq, workflow};
use eframe::egui::{self, Color32};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use theme::color;
use viz::{CutSite, FeatureArc};

fn main() -> eframe::Result<()> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 860.0])
            .with_min_inner_size([1100.0, 720.0])
            .with_title("ClonoDoc"),
        ..Default::default()
    };
    eframe::run_native(
        "clonodoc",
        options,
        Box::new(|cc| {
            theme::apply(&cc.egui_ctx);
            Ok(Box::new(App::new()))
        }),
    )
}

#[derive(PartialEq, Eq, Clone, Copy, Default)]
enum Tab {
    #[default]
    Library,
    Files,
    Check,
    Wetlab,
}

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

/// All per-project (campaign) state. The library is global and lives on `App`;
/// everything here is swapped when the user switches projects in the sidebar.
#[derive(Clone)]
struct ProjectState {
    name: String,
    heavy_vector: Option<String>,
    ground_truth: Vec<GroundTruthRow>,
    order: Vec<SeqRecord>,
    reads: Vec<SeqRecord>,
    order_path: Option<PathBuf>,
    gt_path: Option<PathBuf>,
    gt_is_fasta: bool,
    gt_headers: Vec<String>,
    gt_ab_col: Option<String>,
    gt_heavy_col: Option<String>,
    gt_light_col: Option<String>,
    reads_files: Vec<PathBuf>,
    reads_dir: Option<PathBuf>,
    selected_file: Option<usize>,
    has_overhangs: bool,
    partial_sanger: bool,
    gate1: Vec<Gate1Verdict>,
    gate2: Vec<Gate2Verdict>,
    sel_g1: Option<usize>,
    sel_g2: Option<usize>,
}

impl ProjectState {
    fn new(name: impl Into<String>) -> Self {
        ProjectState {
            name: name.into(),
            heavy_vector: None,
            ground_truth: Vec::new(),
            order: Vec::new(),
            reads: Vec::new(),
            order_path: None,
            gt_path: None,
            gt_is_fasta: false,
            gt_headers: Vec::new(),
            gt_ab_col: None,
            gt_heavy_col: None,
            gt_light_col: None,
            reads_files: Vec::new(),
            reads_dir: None,
            selected_file: None,
            has_overhangs: true,
            partial_sanger: false,
            gate1: Vec::new(),
            gate2: Vec::new(),
            sel_g1: None,
            sel_g2: None,
        }
    }
}

#[derive(Default)]
struct App {
    library: Option<Library>,
    library_path: Option<PathBuf>,
    library_dirty: bool,
    selected_vector: Option<String>,

    // Live working copy of the active project's fields (mirrors ProjectState).
    heavy_vector: Option<String>,
    ground_truth: Vec<GroundTruthRow>,
    order: Vec<SeqRecord>,
    reads: Vec<SeqRecord>,

    order_path: Option<PathBuf>,
    gt_path: Option<PathBuf>,
    gt_is_fasta: bool,
    gt_headers: Vec<String>,
    gt_ab_col: Option<String>,
    gt_heavy_col: Option<String>,
    gt_light_col: Option<String>,
    reads_files: Vec<PathBuf>,
    reads_dir: Option<PathBuf>,
    selected_file: Option<usize>,

    has_overhangs: bool,
    partial_sanger: bool,

    add_vector: AddVectorForm,

    gate1: Vec<Gate1Verdict>,
    gate2: Vec<Gate2Verdict>,
    sel_g1: Option<usize>,
    sel_g2: Option<usize>,

    // Saved project slots; `active_project` indexes the one mirrored above.
    projects: Vec<ProjectState>,
    active_project: usize,

    tab: Tab,
    status: String,
}

impl App {
    fn new() -> Self {
        let mut app = App {
            has_overhangs: true,
            status: "Load or build a library to begin. New here? Library ▸ New, then Add vector from GenBank.".into(),
            projects: vec![ProjectState::new("Project 1")],
            active_project: 0,
            ..Default::default()
        };
        // Dev/demo hook: `CLONODOC_DEMO=<order.fasta>` preloads the bundled library
        // + that order and runs Gate 1, so the rich screens can be exercised
        // headlessly (screenshots / smoke tests). No effect in normal use.
        if let Ok(order_path) = std::env::var("CLONODOC_DEMO") {
            app.load_bundled_library();
            if let Ok(recs) = seqio::read_id_seq_auto(Path::new(&order_path), "") {
                app.order = recs;
                app.order_path = Some(PathBuf::from(&order_path));
            }
            app.run_gate1();
            app.sel_g1 = app.gate1.iter().position(|v| v.passed()).or(Some(0));
            app.tab = Tab::Check;
        }
        // `CLONODOC_DEMO_READS=<reads.fasta>` additionally loads reads, runs Gate 2,
        // and opens the Wetlab tab (for exercising alignment/coverage headlessly).
        if let Ok(reads_path) = std::env::var("CLONODOC_DEMO_READS") {
            if let Ok(recs) = fasta::read_path(Path::new(&reads_path)) {
                app.reads = recs;
                app.reads_files = vec![PathBuf::from(&reads_path)];
            }
            app.run_gate2();
            app.sel_g2 = app.gate2.iter().position(|v| !v.passed()).or(Some(0));
            app.tab = Tab::Wetlab;
        }
        app
    }

    fn project(&self) -> Option<Project> {
        let lib = self.library.as_ref()?;
        Some(workflow::ad_hoc_project(
            lib,
            self.heavy_vector.as_deref(),
            None,
        ))
    }

    // ---- Multi-project state ---------------------------------------------

    /// Copy the live working fields back into the active project slot.
    fn snapshot_active(&mut self) {
        if let Some(slot) = self.projects.get_mut(self.active_project) {
            let name = slot.name.clone();
            *slot = ProjectState {
                name,
                heavy_vector: self.heavy_vector.clone(),
                ground_truth: self.ground_truth.clone(),
                order: self.order.clone(),
                reads: self.reads.clone(),
                order_path: self.order_path.clone(),
                gt_path: self.gt_path.clone(),
                gt_is_fasta: self.gt_is_fasta,
                gt_headers: self.gt_headers.clone(),
                gt_ab_col: self.gt_ab_col.clone(),
                gt_heavy_col: self.gt_heavy_col.clone(),
                gt_light_col: self.gt_light_col.clone(),
                reads_files: self.reads_files.clone(),
                reads_dir: self.reads_dir.clone(),
                selected_file: self.selected_file,
                has_overhangs: self.has_overhangs,
                partial_sanger: self.partial_sanger,
                gate1: self.gate1.clone(),
                gate2: self.gate2.clone(),
                sel_g1: self.sel_g1,
                sel_g2: self.sel_g2,
            };
        }
    }

    /// Load the active project slot into the live working fields.
    fn restore_active(&mut self) {
        let Some(s) = self.projects.get(self.active_project).cloned() else {
            return;
        };
        self.heavy_vector = s.heavy_vector;
        self.ground_truth = s.ground_truth;
        self.order = s.order;
        self.reads = s.reads;
        self.order_path = s.order_path;
        self.gt_path = s.gt_path;
        self.gt_is_fasta = s.gt_is_fasta;
        self.gt_headers = s.gt_headers;
        self.gt_ab_col = s.gt_ab_col;
        self.gt_heavy_col = s.gt_heavy_col;
        self.gt_light_col = s.gt_light_col;
        self.reads_files = s.reads_files;
        self.reads_dir = s.reads_dir;
        self.selected_file = s.selected_file;
        self.has_overhangs = s.has_overhangs;
        self.partial_sanger = s.partial_sanger;
        self.gate1 = s.gate1;
        self.gate2 = s.gate2;
        self.sel_g1 = s.sel_g1;
        self.sel_g2 = s.sel_g2;
    }

    fn switch_project(&mut self, idx: usize) {
        if idx == self.active_project || idx >= self.projects.len() {
            return;
        }
        self.snapshot_active();
        self.active_project = idx;
        self.restore_active();
        self.status = format!("Switched to \u{201c}{}\u{201d}", self.projects[idx].name);
    }

    fn new_project(&mut self) {
        self.snapshot_active();
        let name = format!("Project {}", self.projects.len() + 1);
        self.projects.push(ProjectState::new(name.clone()));
        self.active_project = self.projects.len() - 1;
        self.restore_active();
        // A fresh project defaults its heavy-vector choice to one in the (global) library.
        self.heavy_vector = self
            .library
            .as_ref()
            .and_then(|l| {
                l.vectors
                    .iter()
                    .find(|v| v.chain_class == ChainClass::Heavy)
            })
            .map(|v| v.id.clone());
        self.status = format!("Created \u{201c}{name}\u{201d}. Load its files on Project Files.");
    }

    // ---- Library ---------------------------------------------------------

    fn set_library(&mut self, lib: Library, path: Option<PathBuf>) {
        self.heavy_vector = lib
            .vectors
            .iter()
            .find(|v| v.chain_class == ChainClass::Heavy)
            .map(|v| v.id.clone());
        self.selected_vector = lib.vectors.first().map(|v| v.id.clone());
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

    fn load_bundled_library(&mut self) {
        const BUNDLED: &str = include_str!("../../../reference/example_library.json5");
        match Library::from_json5(BUNDLED) {
            Ok(lib) => self.set_library(lib, None),
            Err(e) => self.status = format!("Bundled library failed to parse: {e}"),
        }
    }

    fn new_library(&mut self) {
        const BUNDLED: &str = include_str!("../../../reference/example_library.json5");
        let mut lib = Library::from_json5(BUNDLED).unwrap_or_else(|_| Library::empty());
        lib.vectors.clear();
        self.set_library(lib, None);
        self.library_dirty = true;
        self.status =
            "New library (overhangs + naming kept, no vectors). Add a vector from GenBank.".into();
    }

    fn save_library_dialog(&mut self) {
        let Some(lib) = &self.library else {
            self.status = "No library to save".into();
            return;
        };
        let mut dialog = rfd::FileDialog::new().add_filter("library", &["json5", "json"]);
        dialog = dialog.set_file_name(
            self.library_path
                .as_ref()
                .and_then(|p| p.file_name().and_then(|n| n.to_str()))
                .unwrap_or("library.json5"),
        );
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
            "Added '{}' ({}, {} bp): insertion {}/{}, anchor {}",
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
        self.selected_vector = Some(v.id.clone());
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
                    "Ground truth: {} antibodies ({} with a sequence)",
                    self.ground_truth.len(),
                    with_seq
                );
                if with_seq == 0 {
                    self.status.push_str(" — check the column mapping.");
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

    fn reload_reads(&mut self) {
        let mut all = Vec::new();
        let mut errors = 0;
        let paths: Vec<PathBuf> = if let Some(dir) = &self.reads_dir {
            list_read_files(dir)
        } else {
            self.reads_files.clone()
        };
        let n_files = paths.len();
        for p in &paths {
            match read_one_reads_file(p) {
                Ok(mut recs) => all.append(&mut recs),
                Err(_) => errors += 1,
            }
        }
        self.reads = all;
        self.status = if errors == 0 {
            format!(
                "Loaded sequencing: {} read(s) from {} file(s)",
                self.reads.len(),
                n_files
            )
        } else {
            format!(
                "Loaded {} read(s) from {} file(s); {errors} failed",
                self.reads.len(),
                n_files
            )
        };
    }

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
    }

    fn run_gate2(&mut self) {
        self.reload_inputs();
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
        // Top spectrum strip (flush to top edge).
        egui::TopBottomPanel::top("spectrum")
            .exact_height(3.0)
            .frame(egui::Frame::new().fill(color::BG_APP))
            .show(ctx, |ui| theme::spectrum_strip(ui, 3.0));

        egui::TopBottomPanel::bottom("status")
            .frame(
                egui::Frame::new()
                    .fill(color::BG_PANEL)
                    .inner_margin(egui::Margin::symmetric(14, 6)),
            )
            .show(ctx, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label(egui::RichText::new("●").color(color::ACCENT).small());
                    ui.label(
                        egui::RichText::new(&self.status)
                            .color(color::TEXT_MUTED)
                            .font(theme::ui_font(12.0)),
                    );
                });
            });

        egui::SidePanel::left("sidebar")
            .exact_width(262.0)
            .frame(
                egui::Frame::new()
                    .fill(color::BG_PANEL)
                    .inner_margin(egui::Margin::same(14)),
            )
            .show(ctx, |ui| self.sidebar(ui));

        egui::CentralPanel::default()
            .frame(egui::Frame::new().fill(color::BG_APP))
            .show(ctx, |ui| {
                self.tab_bar(ui);
                egui::Frame::new()
                    .inner_margin(egui::Margin::same(24))
                    .show(ui, |ui| match self.tab {
                        Tab::Library => self.library_screen(ui),
                        Tab::Files => self.files_screen(ui),
                        Tab::Check => self.check_screen(ui),
                        Tab::Wetlab => self.wetlab_screen(ui),
                    });
            });

        self.add_vector_window(ctx);
    }
}

impl App {
    fn sidebar(&mut self, ui: &mut egui::Ui) {
        // Brand lockup.
        ui.horizontal(|ui| {
            let (rect, _) = ui.allocate_exact_size(egui::vec2(30.0, 18.0), egui::Sense::hover());
            let p = ui.painter();
            p.circle_stroke(
                rect.center() - egui::vec2(6.0, 0.0),
                7.0,
                egui::Stroke::new(2.4, color::ACCENT),
            );
            p.circle_stroke(
                rect.center() + egui::vec2(6.0, 0.0),
                7.0,
                egui::Stroke::new(2.4, color::GOLD),
            );
            ui.vertical(|ui| {
                ui.label(
                    egui::RichText::new("ClonoDoc")
                        .font(theme::ui_font(14.5))
                        .color(color::TEXT_PRIMARY)
                        .strong(),
                );
                theme::section_label(ui, "Cloning Workbench");
            });
        });
        ui.add_space(12.0);

        // Library actions.
        theme::section_label(ui, "Reference Library");
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            if theme::ghost_button(ui, "Load…").clicked() {
                self.load_library_dialog();
            }
            if theme::ghost_button(ui, "New").clicked() {
                self.new_library();
            }
            if theme::ghost_button(ui, "Example").clicked() {
                self.load_bundled_library();
            }
        });
        let has_lib = self.library.is_some();
        ui.add_enabled_ui(has_lib, |ui| {
            if theme::ghost_button(ui, "➕ Add vector from GenBank…").clicked() {
                self.open_add_vector();
            }
            let save = if self.library_dirty {
                "💾 Save library *"
            } else {
                "💾 Save library"
            };
            if theme::ghost_button(ui, save).clicked() {
                self.save_library_dialog();
            }
        });

        if let Some(lib) = &self.library {
            ui.add_space(6.0);
            for v in &lib.vectors {
                let selected = self.selected_vector.as_deref() == Some(v.id.as_str());
                let dot = chain_color(v.chain_class);
                let resp = ui.horizontal(|ui| {
                    theme::type_dot(ui, dot);
                    let txt = egui::RichText::new(&v.display_name)
                        .font(theme::ui_font(13.0))
                        .color(if selected {
                            color::ACCENT
                        } else {
                            color::TEXT_SECONDARY
                        });
                    ui.add(egui::Label::new(txt).sense(egui::Sense::click()))
                        .clicked()
                });
                if resp.inner {
                    self.selected_vector = Some(v.id.clone());
                    self.tab = Tab::Library;
                }
            }
            if lib.vectors.is_empty() {
                ui.label(
                    egui::RichText::new("No vectors yet")
                        .italics()
                        .color(color::TEXT_FAINT)
                        .small(),
                );
            }
        }

        ui.add_space(14.0);
        ui.horizontal(|ui| {
            theme::section_label(ui, "Projects");
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui
                    .add(
                        egui::Label::new(
                            egui::RichText::new("+ New")
                                .font(theme::ui_font(11.0))
                                .color(color::ACCENT),
                        )
                        .sense(egui::Sense::click()),
                    )
                    .on_hover_text(
                        "Create another project (its own ground truth / order / sequencing)",
                    )
                    .clicked()
                {
                    self.new_project();
                }
            });
        });
        ui.add_space(4.0);

        // Project switcher: one selectable row per project, active highlighted.
        let active = self.active_project;
        let names: Vec<String> = self.projects.iter().map(|p| p.name.clone()).collect();
        let mut switch_to = None;
        for (i, name) in names.iter().enumerate() {
            let is_active = i == active;
            let resp = ui.horizontal(|ui| {
                let glyph = if is_active { "▸" } else { "  " };
                let txt = egui::RichText::new(format!("{glyph} {name}"))
                    .font(theme::ui_font(13.0))
                    .color(if is_active {
                        color::ACCENT
                    } else {
                        color::TEXT_SECONDARY
                    });
                ui.add(egui::Label::new(txt).sense(egui::Sense::click()))
                    .clicked()
            });
            if resp.inner && !is_active {
                switch_to = Some(i);
            }
        }
        if let Some(i) = switch_to {
            self.switch_project(i);
        }

        // The active project's loaded inputs.
        ui.add_space(6.0);
        ui.indent("active_inputs", |ui| {
            sidebar_input_row(
                ui,
                "Ground truth",
                self.gt_path.is_some(),
                self.ground_truth.len(),
            );
            sidebar_input_row(ui, "IDT order", self.order_path.is_some(), self.order.len());
            let reads_loaded = self.reads_dir.is_some() || !self.reads_files.is_empty();
            sidebar_input_row(ui, "Sequencing", reads_loaded, self.reads.len());
        });

        ui.add_space(10.0);
        ui.separator();
        ui.checkbox(
            &mut self.has_overhangs,
            egui::RichText::new("Order has overhangs")
                .font(theme::ui_font(12.0))
                .color(color::TEXT_BODY),
        );
        ui.checkbox(
            &mut self.partial_sanger,
            egui::RichText::new("Partial Sanger reads")
                .font(theme::ui_font(12.0))
                .color(color::TEXT_BODY),
        );
    }

    fn tab_bar(&mut self, ui: &mut egui::Ui) {
        egui::Frame::new()
            .fill(color::BG_APP)
            .inner_margin(egui::Margin {
                left: 18,
                right: 18,
                top: 10,
                bottom: 0,
            })
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    self.tab_button(ui, Tab::Library, "Reference Library");
                    self.tab_button(ui, Tab::Files, "Project Files");
                    self.tab_button(ui, Tab::Check, "In-silico Check");
                    self.tab_button(ui, Tab::Wetlab, "Wetlab Verify");
                });
            });
        let (rect, _) =
            ui.allocate_exact_size(egui::vec2(ui.available_width(), 1.0), egui::Sense::hover());
        ui.painter().rect_filled(rect, 0.0, color::BORDER_SOFT);
    }

    fn tab_button(&mut self, ui: &mut egui::Ui, tab: Tab, label: &str) {
        let active = self.tab == tab;
        let color = if active {
            color::ACCENT
        } else {
            Color32::from_rgb(0x8b, 0x98, 0xa5)
        };
        let resp = ui.add(
            egui::Label::new(
                egui::RichText::new(label)
                    .font(theme::ui_font(13.0))
                    .color(color)
                    .strong(),
            )
            .sense(egui::Sense::click()),
        );
        if active {
            let r = resp.rect;
            ui.painter().rect_filled(
                egui::Rect::from_min_max(
                    egui::pos2(r.left(), r.bottom() + 7.0),
                    egui::pos2(r.right(), r.bottom() + 9.0),
                ),
                0.0,
                color::ACCENT,
            );
        }
        if resp.clicked() {
            self.tab = tab;
        }
        ui.add_space(16.0);
    }

    // ---- Screen 1: Reference Library ------------------------------------

    fn library_screen(&mut self, ui: &mut egui::Ui) {
        header(
            ui,
            "Reference Library",
            "Curated backbones & overhangs · reused across projects",
            |ui| {
                if theme::primary_button(ui, "➕ Add vector").clicked() {
                    self.open_add_vector();
                }
                if theme::ghost_button(ui, "Import…").clicked() {
                    self.load_library_dialog();
                }
            },
        );

        let Some(lib) = self.library.clone() else {
            empty_hint(
                ui,
                "No library loaded. Use the sidebar: Load, New, or Example.",
            );
            return;
        };

        detail_panel(ui, "lib_detail", 332.0, |ui| {
            theme::card(ui, |ui| {
                if let Some(v) = self
                    .selected_vector
                    .as_deref()
                    .and_then(|id| lib.vector(id))
                {
                    self.vector_detail(ui, v);
                } else {
                    ui.label(egui::RichText::new("Select a vector").color(color::TEXT_FAINT));
                }
            });
        });
        egui::ScrollArea::vertical()
            .id_salt("lib_list")
            .show(ui, |ui| {
                if lib.vectors.is_empty() {
                    empty_hint(
                        ui,
                        "This library has no vectors yet. Add one from a GenBank file.",
                    );
                }
                for v in &lib.vectors {
                    let selected = self.selected_vector.as_deref() == Some(v.id.as_str());
                    self.library_row(ui, v, selected);
                    ui.add_space(8.0);
                }
            });
    }

    fn library_row(&mut self, ui: &mut egui::Ui, v: &Vector, selected: bool) {
        let frame = egui::Frame::new()
            .fill(if selected {
                color::SEL_BG
            } else {
                color::BG_PANEL
            })
            .stroke(egui::Stroke::new(
                1.0,
                if selected {
                    color::ACCENT.linear_multiply(0.3)
                } else {
                    color::BORDER_CARD
                },
            ))
            .corner_radius(egui::CornerRadius::same(11))
            .inner_margin(egui::Margin::symmetric(16, 14));
        let resp = frame.show(ui, |ui| {
            ui.horizontal(|ui| {
                theme::type_dot(ui, chain_color(v.chain_class));
                ui.add_space(6.0);
                ui.vertical(|ui| {
                    ui.label(
                        egui::RichText::new(&v.display_name)
                            .font(theme::ui_font(14.0))
                            .color(color::TEXT_PRIMARY),
                    );
                    ui.label(
                        egui::RichText::new(format!(
                            "{} · {}",
                            v.chain_class.as_str().to_uppercase(),
                            v.isotype
                        ))
                        .font(theme::ui_font(11.0))
                        .color(color::TEXT_FAINT),
                    );
                });
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        egui::RichText::new(format!("{} bp", viz::group_thousands(v.length)))
                            .font(theme::mono(11.5))
                            .color(color::TEXT_BODY),
                    );
                });
            });
        });
        if resp.response.interact(egui::Sense::click()).clicked() {
            self.selected_vector = Some(v.id.clone());
        }
    }

    fn vector_detail(&self, ui: &mut egui::Ui, v: &Vector) {
        theme::section_label(ui, &format!("{} · selected", v.chain_class.as_str()));
        ui.label(
            egui::RichText::new(&v.display_name)
                .font(theme::ui_font(18.0))
                .color(color::TEXT_PRIMARY)
                .strong(),
        );
        ui.label(
            egui::RichText::new(format!("{} bp", viz::group_thousands(v.length)))
                .font(theme::mono(11.0))
                .color(color::TEXT_MUTED),
        );
        ui.add_space(10.0);

        let arcs = vector_feature_arcs(v);
        viz::circular_plasmid(ui, 240.0, v.length, &v.display_name, &arcs, &[]);
        ui.add_space(10.0);
        // Linear feature track (1 … N).
        viz::linear_feature_bar(ui, v.length, &arcs, 14.0);
        ui.horizontal(|ui| {
            ui.label(
                egui::RichText::new("1")
                    .font(theme::mono(10.0))
                    .color(color::TEXT_FAINT),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    egui::RichText::new(viz::group_thousands(v.length))
                        .font(theme::mono(10.0))
                        .color(color::TEXT_FAINT),
                );
            });
        });
        ui.add_space(10.0);

        for f in &arcs {
            ui.horizontal(|ui| {
                theme::type_dot(ui, f.color);
                ui.add_space(4.0);
                ui.label(
                    egui::RichText::new(&f.name)
                        .font(theme::ui_font(12.5))
                        .color(color::TEXT_BODY),
                );
                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    ui.label(
                        egui::RichText::new(format!("{} bp", f.end.saturating_sub(f.start)))
                            .font(theme::mono(11.0))
                            .color(color::TEXT_MUTED),
                    );
                });
            });
        }
        ui.add_space(8.0);
        let src = v
            .provenance
            .source_file
            .clone()
            .unwrap_or_else(|| "—".into());
        ui.label(
            egui::RichText::new(format!("Source · {src}"))
                .font(theme::ui_font(11.0))
                .color(color::TEXT_FAINT),
        );
    }

    // ---- Screen 2: Project Files ----------------------------------------

    fn files_screen(&mut self, ui: &mut egui::Ui) {
        let proj_name = self
            .projects
            .get(self.active_project)
            .map(|p| p.name.clone())
            .unwrap_or_default();
        header(
            ui,
            "Project Files",
            &format!("Projects / {proj_name} · load this campaign's sequence files"),
            |ui| {
                if theme::primary_button(ui, "Load order…").clicked() {
                    self.load_order_dialog();
                }
                if theme::ghost_button(ui, "+ New project").clicked() {
                    self.new_project();
                }
            },
        );

        detail_panel(ui, "files_preview", 332.0, |ui| {
            theme::card(ui, |ui| {
                theme::section_label(ui, "Preview");
                if let Some(r) = self.selected_file.and_then(|i| self.order.get(i)) {
                    ui.label(
                        egui::RichText::new(&r.id)
                            .font(theme::mono(12.5))
                            .color(color::TEXT_PRIMARY),
                    );
                    let st = seq::detect_type(&r.sequence);
                    ui.label(
                        egui::RichText::new(format!("{:?} · {} bp", st, r.sequence.len()))
                            .font(theme::mono(11.0))
                            .color(color::TEXT_MUTED),
                    );
                    ui.add_space(8.0);
                    ui.label(
                        egui::RichText::new("First bases")
                            .font(theme::ui_font(11.5))
                            .color(color::TEXT_FAINT),
                    );
                    let head: String = r.sequence.chars().take(120).collect();
                    egui::Frame::new()
                        .fill(color::BG_PANEL_ALT)
                        .corner_radius(egui::CornerRadius::same(9))
                        .inner_margin(egui::Margin::same(8))
                        .show(ui, |ui| viz::sequence_row(ui, &head, 11.0, 14.0, None));
                } else {
                    ui.label(
                        egui::RichText::new("Select a record to preview").color(color::TEXT_FAINT),
                    );
                }
            });
        });

        egui::ScrollArea::vertical()
            .id_salt("files_list")
            .show(ui, |ui| {
                // Load buttons row.
                ui.horizontal_wrapped(|ui| {
                    if theme::ghost_button(ui, "📂 Ground truth…").clicked() {
                        self.load_ground_truth_dialog();
                    }
                    if theme::ghost_button(ui, "📂 IDT order…").clicked() {
                        self.load_order_dialog();
                    }
                    if theme::ghost_button(ui, "📂 Sequencing files…").clicked() {
                        self.load_reads_files_dialog();
                    }
                    if theme::ghost_button(ui, "📁 Sequencing folder…").clicked() {
                        self.load_reads_folder_dialog();
                    }
                });

                // Ground-truth column mapping.
                if !self.gt_headers.is_empty() {
                    ui.add_space(10.0);
                    theme::card(ui, |ui| {
                        theme::section_label(ui, "Ground-truth columns");
                        let headers = self.gt_headers.clone();
                        let mut changed = false;
                        changed |= column_picker(
                            ui,
                            "Antibody id",
                            "gt_ab",
                            &headers,
                            &mut self.gt_ab_col,
                        );
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
                    });
                }

                ui.add_space(10.0);
                theme::section_label(ui, &format!("Loaded order records · {}", self.order.len()));
                ui.add_space(4.0);
                let profile = self
                    .library
                    .as_ref()
                    .and_then(|l| l.naming_profiles.first().cloned())
                    .unwrap_or_else(naming::default_profile);
                for (i, r) in self.order.iter().enumerate() {
                    let p = naming::parse_name(&r.id, &profile);
                    let selected = self.selected_file == Some(i);
                    let frame = egui::Frame::new()
                        .fill(if selected {
                            color::SEL_BG
                        } else {
                            color::BG_PANEL
                        })
                        .stroke(egui::Stroke::new(1.0, color::BORDER_CARD))
                        .corner_radius(egui::CornerRadius::same(11))
                        .inner_margin(egui::Margin::symmetric(16, 12));
                    let resp = frame.show(ui, |ui| {
                        ui.horizontal(|ui| {
                            theme::type_dot(ui, chain_color(p.chain_class));
                            ui.add_space(6.0);
                            ui.label(
                                egui::RichText::new(&r.id)
                                    .font(theme::mono(12.5))
                                    .color(color::TEXT_PRIMARY),
                            );
                            ui.with_layout(
                                egui::Layout::right_to_left(egui::Align::Center),
                                |ui| {
                                    role_pill(ui, p.chain_class);
                                    ui.label(
                                        egui::RichText::new(format!("{} bp", r.sequence.len()))
                                            .font(theme::mono(11.0))
                                            .color(color::TEXT_MUTED),
                                    );
                                },
                            );
                        });
                    });
                    if resp.response.interact(egui::Sense::click()).clicked() {
                        self.selected_file = Some(i);
                    }
                    ui.add_space(8.0);
                }
            });
    }

    // ---- Screen 3: In-silico Check (Gate 1) -----------------------------

    fn check_screen(&mut self, ui: &mut egui::Ui) {
        let passed = self.gate1.iter().filter(|v| v.passed()).count();
        let total = self.gate1.len();
        header(
            ui,
            "In-silico Check",
            "Order QC · verify the design before ordering",
            |ui| {
                if theme::primary_button(ui, "▶ Run Gate 1").clicked() {
                    self.run_gate1();
                }
                if !self.gate1.is_empty() && theme::ghost_button(ui, "Export report…").clicked() {
                    self.export_dialog();
                }
                if total > 0 {
                    let c = if passed == total {
                        color::GREEN
                    } else {
                        color::GOLD
                    };
                    theme::pill(ui, &format!("{passed} / {total} passed"), c);
                }
            },
        );

        if self.gate1.is_empty() {
            empty_hint(
                ui,
                "Load a library + IDT order, then Run Gate 1 to verify each ordered construct.",
            );
            return;
        }

        detail_panel(ui, "g1_detail", 400.0, |ui| {
            if let Some(v) = self.sel_g1.and_then(|i| self.gate1.get(i)).cloned() {
                self.gate1_detail(ui, &v);
            } else {
                theme::card(ui, |ui| {
                    ui.label(egui::RichText::new("Select a construct").color(color::TEXT_FAINT));
                    ui.label(
                        egui::RichText::new(
                            "Click a row to see its plasmid map, checks and junction.",
                        )
                        .color(color::TEXT_FAINT)
                        .small(),
                    );
                });
            }
        });
        egui::ScrollArea::vertical()
            .id_salt("g1_list")
            .show(ui, |ui| self.gate1_table(ui));
    }

    fn gate1_table(&mut self, ui: &mut egui::Ui) {
        theme::card(ui, |ui| {
            theme::section_label(ui, "Constructs");
            ui.add_space(6.0);
            egui::Grid::new("g1")
                .striped(true)
                .num_columns(4)
                .spacing(egui::vec2(12.0, 8.0))
                .show(ui, |ui| {
                    for s in ["Record", "Chain", "Verdict", "Read-thru"] {
                        ui.label(
                            egui::RichText::new(s)
                                .font(theme::ui_font(11.0))
                                .color(color::TEXT_FAINT),
                        );
                    }
                    ui.end_row();
                    let mut click = None;
                    for (i, v) in self.gate1.iter().enumerate() {
                        let sel = self.sel_g1 == Some(i);
                        if ui
                            .add(
                                egui::Label::new(
                                    egui::RichText::new(short(&v.record_id, 28))
                                        .font(theme::mono(12.0))
                                        .color(if sel {
                                            color::ACCENT
                                        } else {
                                            color::TEXT_SECONDARY
                                        }),
                                )
                                .sense(egui::Sense::click()),
                            )
                            .clicked()
                        {
                            click = Some(i);
                        }
                        ui.label(
                            egui::RichText::new(&v.chain_class)
                                .font(theme::ui_font(12.0))
                                .color(color::TEXT_BODY),
                        );
                        verdict_chip(ui, v.kind.label());
                        ui.label(
                            egui::RichText::new(opt_bool(v.reads_through))
                                .font(theme::mono(11.5))
                                .color(color::TEXT_MUTED),
                        );
                        ui.end_row();
                    }
                    if click.is_some() {
                        self.sel_g1 = click;
                    }
                });
        });
    }

    fn gate1_detail(&self, ui: &mut egui::Ui, v: &Gate1Verdict) {
        let view = self.build_construct(v);
        theme::card(ui, |ui| {
            theme::section_label(ui, "Simulated construct map");
            ui.label(
                egui::RichText::new(&v.record_id)
                    .font(theme::mono(13.0))
                    .color(color::TEXT_PRIMARY),
            );
            ui.add_space(8.0);
            if let Some(view) = &view {
                viz::circular_plasmid(
                    ui,
                    260.0,
                    view.total_bp,
                    "construct",
                    &view.arcs,
                    &view.cuts,
                );
                ui.add_space(10.0);
                // stat tiles
                ui.horizontal(|ui| {
                    stat_tile(ui, &viz::group_thousands(view.total_bp), "total bp");
                    stat_tile(ui, &view.insert_bp.to_string(), "insert bp");
                    stat_tile(ui, &format!("{:.0}%", view.gc * 100.0), "GC");
                    stat_tile(ui, &view.orf_aa.to_string(), "ORF aa");
                });
            } else {
                ui.label(
                    egui::RichText::new("Map needs a heavy vector + the matching order record.")
                        .color(color::TEXT_FAINT)
                        .small(),
                );
            }
        });
        ui.add_space(12.0);
        theme::card(ui, |ui| {
            theme::section_label(ui, "Verification");
            ui.add_space(4.0);
            check_row(ui, v.passed(), "Verdict", v.kind.label());
            check_row(
                ui,
                v.reads_through == Some(true),
                "Reads through into constant",
                &opt_bool(v.reads_through),
            );
            if let Some(c) = v.core_len {
                check_row(
                    ui,
                    c.is_multiple_of(3),
                    "Insert in frame",
                    &format!("{c} nt"),
                );
            }
            // Only surface a stop as "premature" when the ORF does NOT read through;
            // a stop on a passing construct is the natural constant-region stop.
            if v.reads_through == Some(false) {
                if let Some(s) = v.premature_stop_aa {
                    check_row(ui, false, "Premature stop", &format!("aa {s}"));
                }
            }
            ui.add_space(6.0);
            ui.label(
                egui::RichText::new(&v.reason)
                    .font(theme::ui_font(12.0))
                    .color(color::TEXT_MUTED),
            );
        });
        if let Some(view) = &view {
            if let Some((seq, divider)) = &view.junction {
                ui.add_space(12.0);
                theme::card(ui, |ui| {
                    theme::section_label(ui, "5′ junction — vector ↔ insert");
                    ui.add_space(6.0);
                    base_legend(ui);
                    ui.add_space(4.0);
                    viz::sequence_row(ui, seq, 13.0, 15.0, Some(*divider));
                });
            }
        }
    }

    /// Recompute an assembled-construct view for a Gate-1 row (map + junction + stats).
    fn build_construct(&self, v: &Gate1Verdict) -> Option<ConstructView> {
        let lib = self.library.as_ref()?;
        let project = self.project()?;
        let set = lib.overhang_set(&project.overhang_set)?;
        let rec = self.order.iter().find(|r| r.id == v.record_id)?;
        let class = match v.chain_class.as_str() {
            "heavy" => ChainClass::Heavy,
            "kappa" => ChainClass::Kappa,
            "lambda" => ChainClass::Lambda,
            _ => ChainClass::Light,
        };
        let locus = match class {
            ChainClass::Heavy => Locus::Igh,
            ChainClass::Kappa => Locus::Igk,
            ChainClass::Lambda => Locus::Igl,
            _ => assemble::detect_locus_from_overhangs(
                &rec.sequence,
                set,
                lib.alignment.overhang_max_mismatch,
            )?,
        };
        let vector = project
            .vector_for(locus.chain_class())
            .and_then(|id| lib.vector(id))?;
        let core = if self.has_overhangs {
            let s = assemble::strip_overhangs(
                &rec.sequence,
                set.oh5(locus).unwrap_or(""),
                set.oh3(locus).unwrap_or(""),
                lib.alignment.overhang_max_mismatch,
            );
            if s.oh5_present && s.oh3_present {
                s.core
            } else {
                return None;
            }
        } else {
            seq::clean(&rec.sequence)
        };
        let assembled = assemble::assemble(vector, &core);
        let oh5 = vector.insertion_site.oh5_end;
        let total_bp = assembled.len();
        let insert_bp = core.len();
        let gc = seq::gc_fraction(&core);
        let orf = seq::translate_to_stop(&assembled);
        let orf_aa = orf.len();

        // arcs: insert window (accent) + the vector's own features shifted past the insert.
        let mut arcs = vec![FeatureArc {
            name: "insert".into(),
            start: oh5,
            end: oh5 + insert_bp,
            color: color::ACCENT,
        }];
        for f in vector_feature_arcs(vector) {
            // shift any feature beyond the insertion site by the insert length.
            let shift = |x: usize| {
                if x >= vector.insertion_site.oh3_start {
                    x + insert_bp
                } else {
                    x
                }
            };
            arcs.push(FeatureArc {
                name: f.name,
                start: shift(f.start),
                end: shift(f.end),
                color: f.color,
            });
        }

        // 5' junction window.
        let win = 22usize;
        let start = oh5.saturating_sub(win);
        let end = (oh5 + win).min(assembled.len());
        let junction = if end > start {
            Some((assembled[start..end].to_string(), oh5 - start))
        } else {
            None
        };

        Some(ConstructView {
            total_bp,
            insert_bp,
            gc,
            orf_aa,
            arcs,
            cuts: Vec::new(),
            junction,
        })
    }

    // ---- Screen 4: Wetlab Verify (Gate 2) -------------------------------

    fn wetlab_screen(&mut self, ui: &mut egui::Ui) {
        let pass = self.gate2.iter().filter(|v| v.passed()).count();
        let total = self.gate2.len();
        header(
            ui,
            "Wetlab Verify",
            "Sequencing QC · confirm the finished clone",
            |ui| {
                if theme::primary_button(ui, "▶ Run Gate 2").clicked() {
                    self.run_gate2();
                }
                if !self.gate2.is_empty() && theme::ghost_button(ui, "Export report…").clicked() {
                    self.export_dialog();
                }
                if total > 0 {
                    let c = if pass == total {
                        color::GREEN
                    } else {
                        color::GOLD
                    };
                    theme::pill(ui, &format!("{pass} / {total} pass"), c);
                }
            },
        );

        if self.gate2.is_empty() {
            empty_hint(
                ui,
                "Load sequencing reads (files or a folder), then Run Gate 2.",
            );
            return;
        }

        detail_panel(ui, "g2_detail", 400.0, |ui| {
            if let Some(v) = self.sel_g2.and_then(|i| self.gate2.get(i)).cloned() {
                self.gate2_detail(ui, &v);
            } else {
                theme::card(ui, |ui| {
                    ui.label(egui::RichText::new("Select a read").color(color::TEXT_FAINT));
                });
            }
        });
        egui::ScrollArea::vertical()
            .id_salt("g2_list")
            .show(ui, |ui| self.gate2_table(ui));
    }

    fn gate2_table(&mut self, ui: &mut egui::Ui) {
        theme::card(ui, |ui| {
            theme::section_label(ui, "Reads");
            ui.add_space(6.0);
            egui::Grid::new("g2")
                .striped(true)
                .num_columns(4)
                .spacing(egui::vec2(12.0, 8.0))
                .show(ui, |ui| {
                    for s in ["Record", "Verdict", "Backbone", "Identity"] {
                        ui.label(
                            egui::RichText::new(s)
                                .font(theme::ui_font(11.0))
                                .color(color::TEXT_FAINT),
                        );
                    }
                    ui.end_row();
                    let mut click = None;
                    for (i, v) in self.gate2.iter().enumerate() {
                        let sel = self.sel_g2 == Some(i);
                        if ui
                            .add(
                                egui::Label::new(
                                    egui::RichText::new(short(&v.record_id, 26))
                                        .font(theme::mono(12.0))
                                        .color(if sel {
                                            color::ACCENT
                                        } else {
                                            color::TEXT_SECONDARY
                                        }),
                                )
                                .sense(egui::Sense::click()),
                            )
                            .clicked()
                        {
                            click = Some(i);
                        }
                        verdict_chip(ui, v.kind.label());
                        ui.label(
                            egui::RichText::new(short(
                                v.backbone_vector.as_deref().unwrap_or("—"),
                                14,
                            ))
                            .font(theme::ui_font(11.5))
                            .color(color::TEXT_BODY),
                        );
                        ui.label(
                            egui::RichText::new(
                                v.backbone_identity
                                    .map(|x| format!("{:.1}%", x * 100.0))
                                    .unwrap_or_else(|| "—".into()),
                            )
                            .font(theme::mono(11.5))
                            .color(color::TEXT_MUTED),
                        );
                        ui.end_row();
                    }
                    if click.is_some() {
                        self.sel_g2 = click;
                    }
                });
        });
    }

    fn gate2_detail(&self, ui: &mut egui::Ui, v: &Gate2Verdict) {
        // verdict banner
        let (banner_c, msg) = if v.passed() {
            (color::GREEN, "Construct confirmed")
        } else {
            (color::GOLD, "Needs review")
        };
        egui::Frame::new()
            .fill(banner_c.linear_multiply(0.12))
            .stroke(egui::Stroke::new(1.0, banner_c))
            .corner_radius(egui::CornerRadius::same(14))
            .inner_margin(egui::Margin::same(16))
            .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.label(
                        egui::RichText::new("✓")
                            .font(theme::ui_font(20.0))
                            .color(banner_c),
                    );
                    ui.vertical(|ui| {
                        ui.label(
                            egui::RichText::new(format!("{msg} — {}", v.kind.label()))
                                .font(theme::ui_font(15.0))
                                .color(color::TEXT_PRIMARY),
                        );
                        ui.label(
                            egui::RichText::new(&v.reason)
                                .font(theme::ui_font(12.0))
                                .color(color::TEXT_MUTED),
                        );
                    });
                });
            });
        ui.add_space(12.0);

        theme::card(ui, |ui| {
            theme::section_label(ui, "Result");
            ui.add_space(4.0);
            check_row(
                ui,
                v.backbone_vector.is_some(),
                "Backbone",
                v.backbone_vector.as_deref().unwrap_or("—"),
            );
            check_row(
                ui,
                v.backbone_identity.unwrap_or(0.0) > 0.9,
                "Backbone identity",
                &v.backbone_identity
                    .map(|x| format!("{:.1}%", x * 100.0))
                    .unwrap_or_else(|| "—".into()),
            );
            check_row(
                ui,
                v.reads_through == Some(true),
                "Reads through",
                &opt_bool(v.reads_through),
            );
            check_row(
                ui,
                v.mutations.is_empty(),
                "Mutations",
                &format!("{}", v.mutations.len()),
            );
            check_row(
                ui,
                v.suspected_identity.is_none(),
                "Backbone observed",
                &v.backbone_observed,
            );
        });

        if !v.mutations.is_empty() {
            ui.add_space(12.0);
            theme::card(ui, |ui| {
                theme::section_label(ui, "Mutations in the variable region");
                ui.add_space(4.0);
                for m in &v.mutations {
                    ui.label(
                        egui::RichText::new(format!("{}{}→{}", m.wt, m.position_aa, m.mut_aa))
                            .font(theme::mono(13.0))
                            .color(color::RED),
                    );
                }
            });
        }
        if let Some(s) = &v.suspected_identity {
            ui.add_space(12.0);
            theme::card(ui, |ui| {
                theme::section_label(ui, "Suspected sample swap");
                ui.label(
                    egui::RichText::new(format!("Best match: {s}"))
                        .font(theme::mono(13.0))
                        .color(color::GOLD),
                );
            });
        }

        // Alignment card — observed ORF vs expected protein, windowed around the
        // first mutation (or the V-region start), mismatches highlighted red.
        if let (Some(obs), Some(exp)) = (&v.aligned_observed, &v.aligned_expected) {
            ui.add_space(12.0);
            theme::card(ui, |ui| {
                let leader = v.leader_aa_len.unwrap_or(0);
                // Center the window on the first mutation if present, else the V start.
                let focus = v
                    .mutations
                    .first()
                    .map(|m| leader + m.position_aa)
                    .unwrap_or(leader + 1);
                let half = 13usize;
                let start = focus.saturating_sub(half);
                let end = (focus + half).min(exp.chars().count());
                theme::section_label(ui, "Protein alignment · expected vs observed");
                ui.label(
                    egui::RichText::new(format!("aa {}–{}", start + 1, end))
                        .font(theme::mono(11.0))
                        .color(color::TEXT_MUTED),
                );
                ui.add_space(6.0);
                let slice = |s: &str| s.chars().skip(start).take(end - start).collect::<String>();
                viz::alignment_rows(
                    ui,
                    "Expected",
                    &slice(exp),
                    "Observed",
                    &slice(obs),
                    13.0,
                    14.0,
                );
            });
        }

        // Coverage card — per-base read depth across the expected construct, from
        // all reads of this antibody (most informative for partial-Sanger sets).
        if let Some((depths, total_bp, feats, mean)) = self.coverage_for(v) {
            ui.add_space(12.0);
            theme::card(ui, |ui| {
                theme::section_label(ui, "Read coverage across construct");
                ui.label(
                    egui::RichText::new(format!("mean depth {mean:.0}×"))
                        .font(theme::mono(11.0))
                        .color(color::TEXT_MUTED),
                );
                ui.add_space(6.0);
                viz::coverage_chart(ui, &depths, &feats, total_bp, 120.0);
            });
        }
    }

    /// Compute a coverage profile for the antibody behind a Gate-2 verdict: build
    /// its expected assembled reference, then aggregate every loaded read that
    /// parses to the same `ab_id`. Returns (depths, total_bp, feature arcs, mean).
    fn coverage_for(&self, v: &Gate2Verdict) -> Option<(Vec<f32>, usize, Vec<FeatureArc>, f64)> {
        let lib = self.library.as_ref()?;
        let project = self.project()?;
        let set = lib.overhang_set(&project.overhang_set)?;
        let profile = lib
            .naming_profiles
            .first()
            .cloned()
            .unwrap_or_else(naming::default_profile);

        // Determine the construct's locus/vector + expected core from the order.
        let class = match v.chain_class.as_str() {
            "heavy" => ChainClass::Heavy,
            "kappa" => ChainClass::Kappa,
            "lambda" => ChainClass::Lambda,
            _ => ChainClass::Light,
        };
        let cores = workflow::order_cores(&self.order, lib, &project, self.has_overhangs);
        let key_class = |loc: Locus| loc.chain_class();
        // Find the order core for this antibody (try the parsed class, then light loci).
        let mut core_locus = None;
        for loc in [Locus::Igh, Locus::Igk, Locus::Igl] {
            if cores.contains_key(&(v.ab_id.to_ascii_uppercase(), key_class(loc))) {
                if class == ChainClass::Heavy && loc != Locus::Igh {
                    continue;
                }
                core_locus = Some(loc);
                break;
            }
        }
        let locus = core_locus?;
        let core = cores.get(&(v.ab_id.to_ascii_uppercase(), locus.chain_class()))?;
        let vector = project
            .vector_for(locus.chain_class())
            .and_then(|id| lib.vector(id))?;
        let expected = assemble::assemble(vector, core);
        let _ = set;

        // Gather all reads belonging to this antibody.
        let reads: Vec<String> = self
            .reads
            .iter()
            .filter(|r| {
                naming::parse_name(&r.id, &profile)
                    .ab_id
                    .eq_ignore_ascii_case(&v.ab_id)
            })
            .map(|r| r.sequence.clone())
            .collect();
        if reads.is_empty() {
            return None;
        }
        let cov = clonodoc_core::coverage::coverage_profile(&expected, &reads);
        if cov.is_empty() {
            return None;
        }
        let depths: Vec<f32> = cov.depth.iter().map(|&d| d as f32).collect();
        // Feature arcs (insert + shifted vector features) for the track beneath.
        let oh5 = vector.insertion_site.oh5_end;
        let insert_bp = core.len();
        let mut feats = vec![FeatureArc {
            name: "insert".into(),
            start: oh5,
            end: oh5 + insert_bp,
            color: color::ACCENT,
        }];
        for f in vector_feature_arcs(vector) {
            let shift = |x: usize| {
                if x >= vector.insertion_site.oh3_start {
                    x + insert_bp
                } else {
                    x
                }
            };
            feats.push(FeatureArc {
                name: f.name,
                start: shift(f.start),
                end: shift(f.end),
                color: f.color,
            });
        }
        Some((depths, expected.len(), feats, cov.mean()))
    }
}

// ---- add-vector window -----------------------------------------------------

impl App {
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

        egui::Window::new("Add vector from GenBank").open(&mut open).resizable(true).default_width(460.0).show(ctx, |ui| {
            let f = &mut self.add_vector;
            ui.horizontal(|ui| {
                if ui.button("📂 Choose .gb file…").clicked() { do_pick = true; }
                match &f.gb_path {
                    Some(p) => ui.label(p.file_name().and_then(|n| n.to_str()).unwrap_or("")),
                    None => ui.label(egui::RichText::new("no file chosen").weak()),
                };
            });
            if let Some(gb) = &f.gb {
                ui.label(format!("{} — {} bp, {}", gb.name, gb.length, if gb.circular { "circular" } else { "linear" }));
            }
            ui.separator();
            egui::Grid::new("av_form").num_columns(2).show(ui, |ui| {
                ui.label("Vector id"); ui.text_edit_singleline(&mut f.id); ui.end_row();
                ui.label("Display name"); ui.text_edit_singleline(&mut f.display); ui.end_row();
                ui.label("Isotype"); ui.text_edit_singleline(&mut f.isotype); ui.end_row();
                ui.label("Locus / chain");
                egui::ComboBox::from_id_salt("av_locus").selected_text(locus_label(f.locus)).show_ui(ui, |ui| {
                    ui.selectable_value(&mut f.locus, Locus::Igh, locus_label(Locus::Igh));
                    ui.selectable_value(&mut f.locus, Locus::Igk, locus_label(Locus::Igk));
                    ui.selectable_value(&mut f.locus, Locus::Igl, locus_label(Locus::Igl));
                });
                ui.end_row();
                ui.label("Overhang set");
                egui::ComboBox::from_id_salt("av_set").selected_text(if f.overhang_set.is_empty() { "—".into() } else { f.overhang_set.clone() }).show_ui(ui, |ui| {
                    for id in &overhang_ids { ui.selectable_value(&mut f.overhang_set, id.clone(), id); }
                });
                ui.end_row();
            });
            ui.label(egui::RichText::new("Insertion site & constant anchor are computed automatically (overhangs located in the vector).").weak().small());
            if !f.error.is_empty() {
                ui.colored_label(color::RED, &f.error);
            }
            ui.separator();
            ui.horizontal(|ui| {
                if ui.button("Add vector").clicked() { do_commit = true; }
                if ui.button("Cancel").clicked() { do_cancel = true; }
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
        if self.add_vector.open {
            self.add_vector.open = open;
        }
    }
}

struct ConstructView {
    total_bp: usize,
    insert_bp: usize,
    gc: f64,
    orf_aa: usize,
    arcs: Vec<FeatureArc>,
    cuts: Vec<CutSite>,
    junction: Option<(String, usize)>,
}

// ---- free helpers ----------------------------------------------------------

/// A screen header strip: title + subtitle left, action buttons right, hairline below.
fn header(ui: &mut egui::Ui, title: &str, subtitle: &str, right: impl FnOnce(&mut egui::Ui)) {
    ui.horizontal(|ui| {
        ui.vertical(|ui| {
            theme::page_title(ui, title);
            if !subtitle.is_empty() {
                ui.label(
                    egui::RichText::new(subtitle)
                        .font(theme::ui_font(12.5))
                        .color(color::TEXT_MUTED),
                );
            }
        });
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), right);
    });
    ui.add_space(6.0);
    let (rect, _) =
        ui.allocate_exact_size(egui::vec2(ui.available_width(), 1.0), egui::Sense::hover());
    ui.painter().rect_filled(rect, 0.0, color::DIVIDER_HAIR);
    ui.add_space(16.0);
}

/// A fixed-width right detail column (a `SidePanel`) used by the list+detail
/// screens. Scrolls vertically; never overflows the window. The caller supplies
/// its own card(s) inside.
fn detail_panel(ui: &mut egui::Ui, id: &str, width: f32, add: impl FnOnce(&mut egui::Ui)) {
    egui::SidePanel::right(egui::Id::new(id.to_owned()))
        .resizable(false)
        .exact_width(width)
        .frame(egui::Frame::new().inner_margin(egui::Margin {
            left: 16,
            right: 0,
            top: 0,
            bottom: 0,
        }))
        .show_inside(ui, |ui| {
            egui::ScrollArea::vertical()
                .id_salt(("detail", id))
                .show(ui, add);
        });
}

fn chain_color(c: ChainClass) -> Color32 {
    match c {
        ChainClass::Heavy => color::GOLD,
        ChainClass::Kappa => color::FEATURE_GREEN,
        ChainClass::Lambda => color::PURPLE,
        ChainClass::Light => color::ACCENT,
        ChainClass::Unknown => color::TEXT_MUTED,
    }
}

/// Map a vector's role features to colored arcs (1-based spans → 0-based).
fn vector_feature_arcs(v: &Vector) -> Vec<FeatureArc> {
    let role_color = |role: &str| match role {
        "constant_region" => color::GREEN,
        "resistance" => color::GOLD,
        "origin" => color::PURPLE,
        "promoter" => color::FEATURE_GREEN,
        "signal_peptide" => color::TEXT_MUTED,
        _ => color::TEXT_BODY,
    };
    let mut arcs = Vec::new();
    for (role, f) in &v.features {
        let r = f.range0();
        arcs.push(FeatureArc {
            name: role.replace('_', " "),
            start: r.start,
            end: r.end,
            color: role_color(role),
        });
    }
    arcs
}

fn role_pill(ui: &mut egui::Ui, c: ChainClass) {
    let (label, col) = match c {
        ChainClass::Heavy => ("Heavy", color::GOLD),
        ChainClass::Kappa => ("Kappa", color::FEATURE_GREEN),
        ChainClass::Lambda => ("Lambda", color::PURPLE),
        ChainClass::Light => ("Light", color::ACCENT),
        ChainClass::Unknown => ("?", color::TEXT_MUTED),
    };
    egui::Frame::new()
        .fill(col.linear_multiply(0.14))
        .corner_radius(egui::CornerRadius::same(18))
        .inner_margin(egui::Margin::symmetric(10, 2))
        .show(ui, |ui| {
            ui.label(
                egui::RichText::new(label)
                    .font(theme::ui_font(11.0))
                    .color(col),
            );
        });
}

fn verdict_chip(ui: &mut egui::Ui, label: &str) {
    let c = verdict_color(label);
    egui::Frame::new()
        .fill(c.linear_multiply(0.14))
        .corner_radius(egui::CornerRadius::same(7))
        .inner_margin(egui::Margin::symmetric(8, 2))
        .show(ui, |ui| {
            ui.label(
                egui::RichText::new(label)
                    .font(theme::ui_font(11.0))
                    .color(c),
            );
        });
}

fn check_row(ui: &mut egui::Ui, ok: bool, label: &str, detail: &str) {
    ui.horizontal(|ui| {
        let (sym, c) = if ok {
            ("✓", color::GREEN)
        } else {
            ("!", color::GOLD)
        };
        ui.label(egui::RichText::new(sym).color(c).font(theme::ui_font(13.0)));
        ui.label(
            egui::RichText::new(label)
                .font(theme::ui_font(13.0))
                .color(color::TEXT_BODY),
        );
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            ui.label(
                egui::RichText::new(detail)
                    .font(theme::mono(11.5))
                    .color(color::TEXT_MUTED),
            );
        });
    });
}

fn stat_tile(ui: &mut egui::Ui, value: &str, label: &str) {
    egui::Frame::new()
        .fill(color::BG_APP)
        .corner_radius(egui::CornerRadius::same(9))
        .inner_margin(egui::Margin::same(10))
        .show(ui, |ui| {
            ui.vertical(|ui| {
                ui.label(
                    egui::RichText::new(value)
                        .font(theme::mono(18.0))
                        .color(color::TEXT_PRIMARY),
                );
                ui.label(
                    egui::RichText::new(label)
                        .font(theme::ui_font(10.5))
                        .color(color::TEXT_FAINT),
                );
            });
        });
}

fn base_legend(ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        for b in ['A', 'T', 'C', 'G'] {
            ui.label(
                egui::RichText::new(b.to_string())
                    .font(theme::mono(12.0))
                    .color(theme::base_color(b)),
            );
        }
    });
}

fn sidebar_input_row(ui: &mut egui::Ui, label: &str, loaded: bool, n: usize) {
    ui.horizontal(|ui| {
        let (sym, c) = if loaded {
            ("✔", color::GREEN)
        } else {
            ("—", color::TEXT_FAINT)
        };
        ui.label(egui::RichText::new(sym).color(c).small());
        ui.label(
            egui::RichText::new(label)
                .font(theme::ui_font(12.5))
                .color(color::TEXT_BODY),
        );
        if loaded {
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                ui.label(
                    egui::RichText::new(format!("{n}"))
                        .font(theme::mono(11.0))
                        .color(color::TEXT_MUTED),
                );
            });
        }
    });
}

fn empty_hint(ui: &mut egui::Ui, text: &str) {
    ui.add_space(40.0);
    ui.vertical_centered(|ui| {
        ui.label(
            egui::RichText::new(text)
                .font(theme::ui_font(13.0))
                .color(color::TEXT_FAINT),
        );
    });
}

fn opt_bool(b: Option<bool>) -> String {
    b.map(|x| if x { "yes".into() } else { "no".into() })
        .unwrap_or_else(|| "—".into())
}

fn short(s: &str, n: usize) -> String {
    if s.chars().count() <= n {
        s.to_string()
    } else {
        format!("{}…", s.chars().take(n - 1).collect::<String>())
    }
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

fn column_picker(
    ui: &mut egui::Ui,
    label: &str,
    id: &str,
    headers: &[String],
    selected: &mut Option<String>,
) -> bool {
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label(
            egui::RichText::new(format!("{label}:"))
                .font(theme::ui_font(12.0))
                .color(color::TEXT_BODY),
        );
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

fn verdict_color(label: &str) -> Color32 {
    match label {
        "PASS" => color::GREEN,
        "NO_GROUND_TRUTH" | "SILENT_VARIANT" | "GC_WARNING" | "RARE_CODON_WARNING" => color::GOLD,
        _ => color::RED,
    }
}

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
