//! Data-driven visualizations drawn with egui's painter (`design/README.md`
//! "Assets"): circular plasmid map, linear feature bar, base-colored sequence
//! and junction viewers, read alignment, coverage chart, and chromatogram.
//!
//! These are pure rendering helpers — they take plain data and a target `Ui`.

use crate::theme::{self, color};
use eframe::egui::{self, Color32, CornerRadius, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use std::f32::consts::TAU;

/// One arc/segment of a plasmid or linear map.
#[derive(Clone)]
pub struct FeatureArc {
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub color: Color32,
}

/// A cut-site tick to annotate on the circular map.
#[derive(Clone)]
pub struct CutSite {
    pub name: String,
    pub pos: usize,
}

/// Draw a circular plasmid map: a thin backbone ring, thick colored feature
/// arcs, optional cut-site ticks, and a centered name + bp label.
pub fn circular_plasmid(
    ui: &mut egui::Ui,
    size: f32,
    total_bp: usize,
    name: &str,
    features: &[FeatureArc],
    cuts: &[CutSite],
) {
    let (rect, _) = ui.allocate_exact_size(Vec2::splat(size), Sense::hover());
    let painter = ui.painter_at(rect);
    let center = rect.center();
    let radius = size * 0.40;
    let total = total_bp.max(1) as f32;

    // angle: 12 o'clock = -PI/2, clockwise with increasing bp.
    let angle = |bp: f32| -TAU / 4.0 + (bp / total) * TAU;
    let pt = |bp: f32, r: f32| {
        let a = angle(bp);
        Pos2::new(center.x + r * a.cos(), center.y + r * a.sin())
    };

    // Backbone ring.
    ring(
        &painter,
        center,
        radius,
        2.0,
        color::BORDER_STRONG,
        total_bp,
    );

    // Feature arcs (thick, rounded).
    let arc_r = radius;
    for f in features {
        let (s, e) = (f.start.min(total_bp) as f32, f.end.min(total_bp) as f32);
        let (s, e) = if e >= s { (s, e) } else { (s, e + total) };
        arc(&painter, center, arc_r, angle(s), angle(e), 7.0, f.color);
    }

    // Cut-site ticks + labels.
    for c in cuts {
        let p_in = pt(c.pos as f32, radius - 8.0);
        let p_out = pt(c.pos as f32, radius + 8.0);
        painter.line_segment([p_in, p_out], Stroke::new(1.3, color::TEXT_SECONDARY));
        let lp = pt(c.pos as f32, radius + 20.0);
        painter.text(
            lp,
            egui::Align2::CENTER_CENTER,
            &c.name,
            theme::mono(10.0),
            color::TEXT_MUTED,
        );
    }

    // Center label.
    painter.text(
        center - Vec2::new(0.0, 7.0),
        egui::Align2::CENTER_CENTER,
        name,
        theme::ui_font(13.0),
        color::TEXT_PRIMARY,
    );
    painter.text(
        center + Vec2::new(0.0, 10.0),
        egui::Align2::CENTER_CENTER,
        format!("{} bp", group_thousands(total_bp)),
        theme::mono(11.0),
        color::TEXT_MUTED,
    );
}

fn ring(painter: &egui::Painter, center: Pos2, r: f32, w: f32, c: Color32, _bp: usize) {
    let mut pts = Vec::with_capacity(129);
    for i in 0..=128 {
        let a = (i as f32 / 128.0) * TAU;
        pts.push(Pos2::new(center.x + r * a.cos(), center.y + r * a.sin()));
    }
    painter.add(egui::Shape::line(pts, Stroke::new(w, c)));
}

fn arc(painter: &egui::Painter, center: Pos2, r: f32, a0: f32, a1: f32, w: f32, c: Color32) {
    let steps = ((a1 - a0).abs() / TAU * 128.0).ceil().max(2.0) as usize;
    let mut pts = Vec::with_capacity(steps + 1);
    for i in 0..=steps {
        let a = a0 + (a1 - a0) * (i as f32 / steps as f32);
        pts.push(Pos2::new(center.x + r * a.cos(), center.y + r * a.sin()));
    }
    painter.add(egui::Shape::line(pts, Stroke::new(w, c)));
}

/// A horizontal feature bar (linear map): colored segments proportional to span,
/// with optional `1 … N` end labels.
pub fn linear_feature_bar(
    ui: &mut egui::Ui,
    total_bp: usize,
    features: &[FeatureArc],
    height: f32,
) {
    let width = ui.available_width();
    let (rect, _) = ui.allocate_exact_size(Vec2::new(width, height), Sense::hover());
    let painter = ui.painter_at(rect);
    let total = total_bp.max(1) as f32;
    // base track
    painter.rect_filled(rect, CornerRadius::same(4), color::SURFACE_CHIP);
    for f in features {
        let x0 = rect.left() + (f.start.min(total_bp) as f32 / total) * rect.width();
        let x1 = rect.left() + (f.end.min(total_bp) as f32 / total) * rect.width();
        let seg = Rect::from_min_max(
            Pos2::new(x0, rect.top()),
            Pos2::new(x1.max(x0 + 1.0), rect.bottom()),
        );
        painter.rect_filled(seg, CornerRadius::same(3), f.color);
    }
}

/// A base-colored monospace sequence row inside a horizontally scrollable block.
/// `cell` is the per-character advance; `divider_at` optionally draws a dashed
/// accent junction divider before that 0-based index.
pub fn sequence_row(ui: &mut egui::Ui, seq: &str, cell: f32, size: f32, divider_at: Option<usize>) {
    let chars: Vec<char> = seq.chars().collect();
    let width = chars.len() as f32 * cell;
    let height = size + 10.0;
    egui::ScrollArea::horizontal()
        .id_salt(("seqrow", seq.len(), divider_at))
        .show(ui, |ui| {
            let (rect, _) = ui.allocate_exact_size(
                Vec2::new(width.max(ui.available_width()), height),
                Sense::hover(),
            );
            let painter = ui.painter_at(rect);
            for (i, ch) in chars.iter().enumerate() {
                let x = rect.left() + i as f32 * cell + cell / 2.0;
                painter.text(
                    Pos2::new(x, rect.center().y),
                    egui::Align2::CENTER_CENTER,
                    ch.to_string(),
                    FontId::new(size, egui::FontFamily::Name(theme::MONO.into())),
                    theme::base_color(*ch),
                );
            }
            if let Some(d) = divider_at {
                let x = rect.left() + d as f32 * cell;
                dashed_vline(&painter, x, rect.top(), rect.bottom(), color::ACCENT);
            }
        });
}

/// Two stacked sequence rows (expected ref over read) with mismatches boxed red.
///
/// Scaffolding for the Wetlab alignment card — wired once `gate2` surfaces the
/// aligned `expected`/`observed` strings (a small `clonodoc-core` addition).
#[allow(dead_code)]
pub fn alignment_rows(
    ui: &mut egui::Ui,
    ref_label: &str,
    reference: &str,
    read_label: &str,
    read: &str,
    cell: f32,
    size: f32,
) {
    let r: Vec<char> = reference.chars().collect();
    let q: Vec<char> = read.chars().collect();
    let n = r.len().max(q.len());
    let label_w = 120.0;
    let width = label_w + n as f32 * cell;
    egui::ScrollArea::horizontal()
        .id_salt(("align", reference.len(), read.len()))
        .show(ui, |ui| {
            let (rect, _) = ui.allocate_exact_size(
                Vec2::new(width.max(ui.available_width()), size * 2.0 + 24.0),
                Sense::hover(),
            );
            let painter = ui.painter_at(rect);
            let row_y = |k: usize| rect.top() + 8.0 + size / 2.0 + k as f32 * (size + 8.0);
            painter.text(
                Pos2::new(rect.left(), row_y(0)),
                egui::Align2::LEFT_CENTER,
                ref_label,
                theme::ui_font(11.5),
                color::TEXT_MUTED,
            );
            painter.text(
                Pos2::new(rect.left(), row_y(1)),
                egui::Align2::LEFT_CENTER,
                read_label,
                theme::ui_font(11.5),
                color::TEXT_MUTED,
            );
            for i in 0..n {
                let x = rect.left() + label_w + i as f32 * cell + cell / 2.0;
                let rc = r.get(i).copied();
                let qc = q.get(i).copied();
                if let Some(c) = rc {
                    painter.text(
                        Pos2::new(x, row_y(0)),
                        egui::Align2::CENTER_CENTER,
                        c.to_string(),
                        theme::mono(size),
                        theme::base_color(c),
                    );
                }
                if let Some(c) = qc {
                    let mismatch = rc.is_some() && rc != qc;
                    if mismatch {
                        let cellrect = Rect::from_center_size(
                            Pos2::new(x, row_y(1)),
                            Vec2::new(cell - 1.0, size + 4.0),
                        );
                        painter.rect_filled(
                            cellrect,
                            CornerRadius::same(3),
                            color::RED.linear_multiply(0.22),
                        );
                    }
                    let col = if mismatch {
                        color::RED
                    } else {
                        theme::base_color(c)
                    };
                    painter.text(
                        Pos2::new(x, row_y(1)),
                        egui::Align2::CENTER_CENTER,
                        c.to_string(),
                        theme::mono(size),
                        col,
                    );
                }
            }
        });
}

/// A filled coverage-depth area chart with an optional feature bar beneath.
///
/// Scaffolding for the Wetlab coverage card — wired once `gate2` retains
/// per-base read depth.
#[allow(dead_code)]
pub fn coverage_chart(
    ui: &mut egui::Ui,
    depths: &[f32],
    features: &[FeatureArc],
    total_bp: usize,
    height: f32,
) {
    let width = ui.available_width();
    let (rect, _) = ui.allocate_exact_size(Vec2::new(width, height), Sense::hover());
    let painter = ui.painter_at(rect);
    let bar_h = 16.0;
    let chart = Rect::from_min_max(
        rect.min,
        Pos2::new(rect.right(), rect.bottom() - bar_h - 6.0),
    );
    painter.rect_filled(chart, CornerRadius::same(6), color::BG_APP);

    if !depths.is_empty() {
        let max = depths.iter().cloned().fold(1.0_f32, f32::max);
        let n = depths.len();
        let mut poly: Vec<Pos2> = Vec::with_capacity(n + 2);
        poly.push(Pos2::new(chart.left(), chart.bottom()));
        for (i, d) in depths.iter().enumerate() {
            let x = chart.left() + (i as f32 / (n - 1).max(1) as f32) * chart.width();
            let y = chart.bottom() - (d / max) * (chart.height() - 4.0);
            poly.push(Pos2::new(x, y));
        }
        poly.push(Pos2::new(chart.right(), chart.bottom()));
        painter.add(egui::Shape::convex_polygon(
            poly.clone(),
            color::ACCENT.linear_multiply(0.25),
            Stroke::NONE,
        ));
        // top outline
        painter.add(egui::Shape::line(
            poly[1..poly.len() - 1].to_vec(),
            Stroke::new(1.5, color::ACCENT),
        ));
    }

    // feature bar beneath
    let total = total_bp.max(1) as f32;
    let bar = Rect::from_min_max(
        Pos2::new(rect.left(), rect.bottom() - bar_h),
        Pos2::new(rect.right(), rect.bottom()),
    );
    painter.rect_filled(bar, CornerRadius::same(4), color::SURFACE_CHIP);
    for f in features {
        let x0 = bar.left() + (f.start.min(total_bp) as f32 / total) * bar.width();
        let x1 = bar.left() + (f.end.min(total_bp) as f32 / total) * bar.width();
        let seg = Rect::from_min_max(
            Pos2::new(x0, bar.top()),
            Pos2::new(x1.max(x0 + 1.0), bar.bottom()),
        );
        painter.rect_filled(seg, CornerRadius::same(3), f.color);
    }
}

/// A synthetic chromatogram: one bell curve per base (quadratic bezier), colored
/// by base, with a base-letter row beneath and an optional highlighted position.
///
/// Scaffolding for the Wetlab chromatogram card — wired once AB1 traces are
/// retained (v1 intentionally parses base calls only, not trace peaks).
#[allow(dead_code)]
pub fn chromatogram(ui: &mut egui::Ui, bases: &str, highlight: Option<usize>, height: f32) {
    let chars: Vec<char> = bases.chars().collect();
    let cell = 30.0;
    let width = (chars.len() as f32 * cell).max(ui.available_width());
    egui::ScrollArea::horizontal()
        .id_salt(("chroma", bases.len()))
        .show(ui, |ui| {
            let (rect, _) = ui.allocate_exact_size(Vec2::new(width, height), Sense::hover());
            let painter = ui.painter_at(rect);
            let baseline = rect.bottom() - 22.0;
            painter.line_segment(
                [
                    Pos2::new(rect.left(), baseline),
                    Pos2::new(rect.right(), baseline),
                ],
                Stroke::new(1.0, color::BORDER_CARD),
            );
            for (i, ch) in chars.iter().enumerate() {
                let cx = rect.left() + i as f32 * cell + cell / 2.0;
                if highlight == Some(i) {
                    let hl = Rect::from_min_max(
                        Pos2::new(cx - cell / 2.0, rect.top()),
                        Pos2::new(cx + cell / 2.0, baseline),
                    );
                    painter.rect_filled(
                        hl,
                        CornerRadius::same(4),
                        color::GOLD.linear_multiply(0.12),
                    );
                }
                let col = theme::base_color(*ch);
                // bell: quadratic from (cx-14, baseline) control (cx, baseline-peak) to (cx+14, baseline)
                let peak = height - 40.0;
                let p0 = Pos2::new(cx - 14.0, baseline);
                let ctrl = Pos2::new(cx, baseline - peak);
                let p2 = Pos2::new(cx + 14.0, baseline);
                let mut pts = Vec::with_capacity(17);
                for s in 0..=16 {
                    let t = s as f32 / 16.0;
                    let mt = 1.0 - t;
                    let p = Pos2::new(
                        mt * mt * p0.x + 2.0 * mt * t * ctrl.x + t * t * p2.x,
                        mt * mt * p0.y + 2.0 * mt * t * ctrl.y + t * t * p2.y,
                    );
                    pts.push(p);
                }
                painter.add(egui::Shape::line(pts, Stroke::new(2.0, col)));
                painter.text(
                    Pos2::new(cx, baseline + 11.0),
                    egui::Align2::CENTER_CENTER,
                    ch.to_string(),
                    theme::mono(13.0),
                    col,
                );
            }
        });
}

/// Format an integer with thousands separators (e.g. 2686 -> "2,686").
pub fn group_thousands(n: usize) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::new();
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

fn dashed_vline(painter: &egui::Painter, x: f32, top: f32, bottom: f32, c: Color32) {
    let dash = 4.0;
    let gap = 3.0;
    let mut y = top;
    while y < bottom {
        let y1 = (y + dash).min(bottom);
        painter.line_segment([Pos2::new(x, y), Pos2::new(x, y1)], Stroke::new(2.0, c));
        y = y1 + gap;
    }
}
