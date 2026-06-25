//! Visual theme for ClonoDoc — the "Bench" dark design (`design/README.md`).
//!
//! Centralizes the design tokens (colors, fonts, radii, spacing) and a set of
//! reusable widget helpers (cards, chips, pills, buttons, section labels, the
//! top spectrum strip, DNA base colors) so the screens read declaratively and
//! stay consistent. Approximates the hi-fi web mock within egui's capabilities.

use eframe::egui::{self, Color32, CornerRadius, FontFamily, FontId, Stroke};

// ---- Color tokens (design/README.md "Design Tokens") -----------------------

pub mod color {
    use eframe::egui::Color32;
    const fn hex(r: u8, g: u8, b: u8) -> Color32 {
        Color32::from_rgb(r, g, b)
    }
    // surfaces & text
    pub const BG_APP: Color32 = hex(0x0E, 0x12, 0x17);
    pub const BG_PANEL: Color32 = hex(0x11, 0x16, 0x1C);
    pub const BG_PANEL_ALT: Color32 = hex(0x10, 0x15, 0x1B);
    pub const SURFACE_CHIP: Color32 = hex(0x16, 0x1D, 0x24);
    pub const SURFACE_MUTED: Color32 = hex(0x1A, 0x22, 0x2B);
    pub const BORDER_STRONG: Color32 = hex(0x2A, 0x33, 0x3E);
    pub const BORDER_CARD: Color32 = hex(0x1F, 0x27, 0x2F);
    pub const BORDER_SOFT: Color32 = hex(0x23, 0x2B, 0x34);
    pub const DIVIDER_HAIR: Color32 = hex(0x16, 0x1D, 0x24);
    pub const TEXT_PRIMARY: Color32 = hex(0xE7, 0xED, 0xF3);
    pub const TEXT_SECONDARY: Color32 = hex(0xCD, 0xD6, 0xDF);
    pub const TEXT_BODY: Color32 = hex(0x9A, 0xA7, 0xB4);
    pub const TEXT_MUTED: Color32 = hex(0x7E, 0x8C, 0x99);
    pub const TEXT_FAINT: Color32 = hex(0x5E, 0x6B, 0x78);

    // brand & semantic (Scripps palette)
    pub const ACCENT: Color32 = hex(0x5B, 0xC2, 0xD2); // Nitrogen Sky
    pub const GOLD: Color32 = hex(0xFF, 0xC9, 0x51); // Infinity Gold
    pub const GREEN: Color32 = hex(0x57, 0xC9, 0x8A); // pass
    pub const PURPLE: Color32 = hex(0x82, 0x91, 0xC6); // Mobius Blue / ori
    pub const FEATURE_GREEN: Color32 = hex(0x6E, 0xB7, 0x44); // lacZα / enzymes
    pub const RED: Color32 = hex(0xEF, 0x6A, 0x53); // mismatch / error

    // selection tint used on active rows / tabs
    pub const SEL_BG: Color32 = hex(0x13, 0x20, 0x2A);

    // DNA base coding
    pub const BASE_A: Color32 = FEATURE_GREEN;
    pub const BASE_T: Color32 = RED;
    pub const BASE_C: Color32 = ACCENT;
    pub const BASE_G: Color32 = GOLD;
    pub const BASE_N: Color32 = TEXT_BODY;

    /// The 5-stop brand spectrum used for the top strip.
    pub const SPECTRUM: [Color32; 5] = [
        hex(0xFF, 0xC9, 0x51),
        hex(0xEF, 0x63, 0x4F),
        hex(0xA5, 0x44, 0x72),
        hex(0x19, 0x3D, 0x66),
        hex(0x5B, 0xC2, 0xD2),
    ];
}

/// Color for a single DNA base character.
pub fn base_color(b: char) -> Color32 {
    match b.to_ascii_uppercase() {
        'A' => color::BASE_A,
        'T' | 'U' => color::BASE_T,
        'C' => color::BASE_C,
        'G' => color::BASE_G,
        _ => color::BASE_N,
    }
}

// ---- Font families ---------------------------------------------------------

/// Logical mono family name used for sequence text, lengths, positions.
pub const MONO: &str = "plex_mono";

/// Install IBM Plex Sans (proportional) + IBM Plex Mono (monospace) and tune the
/// global text styles to the design's type scale.
pub fn install_fonts(ctx: &egui::Context) {
    use egui::{FontData, FontDefinitions};
    let mut fonts = FontDefinitions::default();

    fonts.font_data.insert(
        "plex_sans".to_owned(),
        FontData::from_static(include_bytes!("../assets/fonts/IBMPlexSans.ttf")).into(),
    );
    fonts.font_data.insert(
        "plex_mono".to_owned(),
        FontData::from_static(include_bytes!("../assets/fonts/IBMPlexMono-Regular.ttf")).into(),
    );

    // IBM Plex Sans becomes the primary proportional face.
    fonts
        .families
        .entry(FontFamily::Proportional)
        .or_default()
        .insert(0, "plex_sans".to_owned());
    // IBM Plex Mono leads the monospace family (used for all sequence text).
    fonts
        .families
        .entry(FontFamily::Monospace)
        .or_default()
        .insert(0, "plex_mono".to_owned());
    // A named family alias so sequence widgets can request mono explicitly.
    fonts
        .families
        .insert(FontFamily::Name(MONO.into()), vec!["plex_mono".to_owned()]);

    ctx.set_fonts(fonts);
}

/// A monospace [`FontId`] at the given size (for sequence/numeric text).
pub fn mono(size: f32) -> FontId {
    FontId::new(size, FontFamily::Name(MONO.into()))
}

/// A proportional [`FontId`] at the given size.
pub fn ui_font(size: f32) -> FontId {
    FontId::new(size, FontFamily::Proportional)
}

// ---- Global visuals --------------------------------------------------------

/// Apply the dark palette, rounding, and spacing to the whole context.
pub fn apply(ctx: &egui::Context) {
    install_fonts(ctx);

    let mut style = (*ctx.style()).clone();
    let mut v = egui::Visuals::dark();

    v.dark_mode = true;
    v.override_text_color = Some(color::TEXT_BODY);
    v.panel_fill = color::BG_APP;
    v.window_fill = color::BG_PANEL;
    v.extreme_bg_color = color::BG_APP; // text-edit / inset backgrounds
    v.faint_bg_color = color::BG_PANEL;
    v.window_stroke = Stroke::new(1.0, color::BORDER_CARD);
    v.selection.bg_fill = color::ACCENT.linear_multiply(0.30);
    v.selection.stroke = Stroke::new(1.0, color::ACCENT);
    v.hyperlink_color = color::ACCENT;

    let radius = CornerRadius::same(9);
    for w in [
        &mut v.widgets.noninteractive,
        &mut v.widgets.inactive,
        &mut v.widgets.hovered,
        &mut v.widgets.active,
        &mut v.widgets.open,
    ] {
        w.corner_radius = radius;
        w.bg_stroke = Stroke::new(1.0, color::BORDER_CARD);
        w.fg_stroke = Stroke::new(1.0, color::TEXT_BODY);
    }
    v.widgets.noninteractive.bg_fill = color::BG_PANEL;
    v.widgets.inactive.bg_fill = color::SURFACE_CHIP;
    v.widgets.inactive.weak_bg_fill = color::SURFACE_CHIP;
    v.widgets.hovered.bg_fill = color::SURFACE_MUTED;
    v.widgets.hovered.weak_bg_fill = color::SURFACE_MUTED;
    v.widgets.hovered.bg_stroke = Stroke::new(1.0, color::ACCENT);
    v.widgets.hovered.fg_stroke = Stroke::new(1.0, color::TEXT_PRIMARY);
    v.widgets.active.bg_fill = color::SURFACE_MUTED;
    v.widgets.active.fg_stroke = Stroke::new(1.0, color::TEXT_PRIMARY);

    style.visuals = v;

    // Type scale (design tokens). egui maps these to text styles.
    use egui::TextStyle::*;
    style.text_styles = [
        (Heading, ui_font(21.0)),
        (Body, ui_font(13.0)),
        (Button, ui_font(13.0)),
        (Small, ui_font(11.0)),
        (Monospace, mono(13.0)),
    ]
    .into();

    style.spacing.item_spacing = egui::vec2(8.0, 8.0);
    style.spacing.button_padding = egui::vec2(12.0, 7.0);
    style.spacing.window_margin = egui::Margin::same(0);

    ctx.set_style(style);
}

// ---- Reusable widgets ------------------------------------------------------

/// Paint the 3px brand spectrum strip across a full-width rect.
pub fn spectrum_strip(ui: &mut egui::Ui, height: f32) {
    let (rect, _) = ui.allocate_exact_size(
        egui::vec2(ui.available_width(), height),
        egui::Sense::hover(),
    );
    let stops = &color::SPECTRUM;
    let n = stops.len();
    let seg = rect.width() / (n - 1) as f32;
    let painter = ui.painter();
    for i in 0..n - 1 {
        let x0 = rect.left() + seg * i as f32;
        let x1 = rect.left() + seg * (i + 1) as f32;
        // Approximate the gradient with a few vertical slices between stops.
        let slices = 24;
        for s in 0..slices {
            let t0 = s as f32 / slices as f32;
            let c = lerp_color(stops[i], stops[i + 1], t0);
            let sx0 = x0 + (x1 - x0) * t0;
            let sx1 = x0 + (x1 - x0) * (s + 1) as f32 / slices as f32;
            painter.rect_filled(
                egui::Rect::from_min_max(
                    egui::pos2(sx0, rect.top()),
                    egui::pos2(sx1, rect.bottom()),
                ),
                0.0,
                c,
            );
        }
    }
}

pub fn lerp_color(a: Color32, b: Color32, t: f32) -> Color32 {
    let l = |x: u8, y: u8| (x as f32 + (y as f32 - x as f32) * t) as u8;
    Color32::from_rgb(l(a.r(), b.r()), l(a.g(), b.g()), l(a.b(), b.b()))
}

/// A card frame: panel fill, 1px card border, 14px radius, generous padding.
pub fn card(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui)) {
    egui::Frame::new()
        .fill(color::BG_PANEL)
        .stroke(Stroke::new(1.0, color::BORDER_CARD))
        .corner_radius(CornerRadius::same(14))
        .inner_margin(egui::Margin::same(18))
        .show(ui, add);
}

/// An uppercase, letter-spaced section label (e.g. "REFERENCE LIBRARY").
pub fn section_label(ui: &mut egui::Ui, text: &str) {
    // egui has no letter-spacing; approximate by inserting thin spaces.
    let spaced: String = text
        .to_uppercase()
        .chars()
        .flat_map(|c| [c, '\u{2009}'])
        .collect();
    ui.label(
        egui::RichText::new(spaced)
            .font(ui_font(10.5))
            .color(color::TEXT_FAINT),
    );
}

/// A page title (21px primary).
pub fn page_title(ui: &mut egui::Ui, text: &str) {
    ui.label(
        egui::RichText::new(text)
            .font(ui_font(21.0))
            .color(color::TEXT_PRIMARY),
    );
}

/// A primary (solid accent, dark text) button.
pub fn primary_button(ui: &mut egui::Ui, text: &str) -> egui::Response {
    let btn = egui::Button::new(
        egui::RichText::new(text)
            .color(color::BG_APP)
            .font(ui_font(13.0)),
    )
    .fill(color::ACCENT)
    .corner_radius(CornerRadius::same(9));
    ui.add(btn)
}

/// A ghost (transparent, bordered) button.
pub fn ghost_button(ui: &mut egui::Ui, text: &str) -> egui::Response {
    let btn = egui::Button::new(
        egui::RichText::new(text)
            .color(color::TEXT_SECONDARY)
            .font(ui_font(13.0)),
    )
    .fill(Color32::TRANSPARENT)
    .stroke(Stroke::new(1.0, color::BORDER_STRONG))
    .corner_radius(CornerRadius::same(9));
    ui.add(btn)
}

/// A status pill (rounded, tinted) — color sets text + a translucent fill.
pub fn pill(ui: &mut egui::Ui, text: &str, c: Color32) {
    egui::Frame::new()
        .fill(c.linear_multiply(0.14))
        .corner_radius(CornerRadius::same(20))
        .inner_margin(egui::Margin::symmetric(12, 4))
        .show(ui, |ui| {
            ui.label(egui::RichText::new(text).color(c).font(ui_font(12.0)));
        });
}

/// A small colored type-dot (10px square, rounded).
pub fn type_dot(ui: &mut egui::Ui, c: Color32) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(10.0, 10.0), egui::Sense::hover());
    ui.painter().rect_filled(rect, CornerRadius::same(3), c);
}
