use std::sync::{Arc, mpsc::Sender};

use aviutl2_eframe::{AviUtl2EframeHandle, eframe, egui};

use crate::{
    model::{MAX_PINS, MIN_PINS, PinKind, PinModel, PinOperation, Point},
    state::{CapturedFrame, SharedState},
    worker::WorkerCommand,
};

const IMAGE_MARGIN: f32 = 16.0;
const PIN_RADIUS: f32 = 6.0;
const PIN_HOVER_RADIUS: f32 = 12.0;
const BONE_LINE_HOVER_RADIUS: f32 = 10.0;
const BONE_ARROW_SIZE: f32 = 7.0;
const MIN_ZOOM: f32 = 0.1;
const MAX_ZOOM: f32 = 16.0;
const HELP_TOOLTIP_DELAY_SECONDS: f64 = 1.0;
const HELP_TEXT: &str = "操作方法\n\
・左クリック: ピンを追加\n\
・ピンを左クリック: ピンを選択\n\
・ピンを左ドラッグ: ピンを移動\n\
・ピンをダブルクリック: 制御フィルタを追加\n\
・Shift+左クリック: 選択ボーンの子として接続\n\
・ボーンの線を左クリック: ボーンピンを挿入\n\
・ピンを右クリック: ピンを削除\n\
・空白を右クリック: 選択解除\n\
・ホイール／Shift+ホイール: 縦／横移動\n\
・Ctrl+ホイール: 拡大縮小\n\
・中ボタンドラッグ: 表示位置を移動";

#[derive(Clone, Copy, Debug)]
struct ViewState {
    zoom: f32,
    pan: egui::Vec2,
}

impl Default for ViewState {
    fn default() -> Self {
        Self {
            zoom: 1.0,
            pan: egui::Vec2::ZERO,
        }
    }
}

pub struct PuppetEditorApp {
    _handle: AviUtl2EframeHandle,
    shared: Arc<SharedState>,
    sender: Sender<WorkerCommand>,
    texture: Option<egui::TextureHandle>,
    texture_sequence: u64,
    model_retry_capture_after: u64,
    local_model: Option<PinModel>,
    model_revision: u64,
    dragging_pin: Option<usize>,
    pin_kind: PinKind,
    selected_pin: Option<usize>,
    last_help_pointer: Option<egui::Pos2>,
    help_hover_started_at: Option<f64>,
    view: ViewState,
}

impl PuppetEditorApp {
    pub fn new(
        cc: &eframe::CreationContext<'_>,
        handle: AviUtl2EframeHandle,
        shared: Arc<SharedState>,
        sender: Sender<WorkerCommand>,
    ) -> Self {
        cc.egui_ctx.all_styles_mut(|style| {
            style.visuals = aviutl2_eframe::aviutl2_visuals();
        });
        cc.egui_ctx.set_fonts(aviutl2_eframe::aviutl2_fonts());
        shared.set_egui_context(cc.egui_ctx.clone());
        Self {
            _handle: handle,
            shared,
            sender,
            texture: None,
            texture_sequence: 0,
            model_retry_capture_after: u64::MAX,
            local_model: None,
            model_revision: 0,
            dragging_pin: None,
            pin_kind: PinKind::Position,
            selected_pin: None,
            last_help_pointer: None,
            help_hover_started_at: None,
            view: ViewState::default(),
        }
    }

    fn refresh_local_state(&mut self) {
        let editor = self.shared.editor_snapshot();
        if editor.model.is_none() {
            self.selected_pin = None;
        }
        if self.dragging_pin.is_none() && !editor.busy && editor.revision != self.model_revision {
            self.local_model = editor.model.clone();
            self.model_revision = editor.revision;
            if self.selected_pin.is_some_and(|selected| {
                self.local_model
                    .as_ref()
                    .is_none_or(|model| selected >= model.src.len())
            }) {
                self.selected_pin = None;
            }
        }
    }

    fn refresh_texture(&mut self, context: &egui::Context) -> Option<CapturedFrame> {
        let editor = self.shared.editor_snapshot();
        let capture = self.shared.capture_snapshot();
        let frame = capture
            .latest
            .filter(|frame| frame.sequence > editor.capture_after)?;
        if frame.sequence != self.texture_sequence {
            let image = egui::ColorImage::from_rgba_unmultiplied(
                [frame.width, frame.height],
                frame.rgba.as_ref(),
            );
            if let Some(texture) = &mut self.texture {
                texture.set(image, egui::TextureOptions::LINEAR);
            } else {
                self.texture = Some(context.load_texture(
                    "puppet-tool-input",
                    image,
                    egui::TextureOptions::LINEAR,
                ));
            }
            self.texture_sequence = frame.sequence;
        }
        if editor.model.is_none() && editor.capture_after != self.model_retry_capture_after {
            self.model_retry_capture_after = editor.capture_after;
            let _ = self.sender.send(WorkerCommand::Refresh);
        }
        Some(frame)
    }

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        let editor = self.shared.editor_snapshot();
        egui::Panel::top("puppet_toolbar").show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                if ui
                    .add_enabled(!editor.busy, egui::Button::new("＋"))
                    .on_hover_text("オブジェクト追加")
                    .clicked()
                {
                    let _ = self.sender.send(WorkerCommand::AddObject);
                }
                let mut effect_number = editor
                    .target
                    .as_ref()
                    .map_or(1, |target| target.effect_index + 1);
                ui.label("効果番号");
                if ui
                    .add(
                        egui::DragValue::new(&mut effect_number)
                            .range(1..=20)
                            .speed(0.1),
                    )
                    .changed()
                {
                    let _ = self.sender.send(WorkerCommand::SelectFocused {
                        effect_index: effect_number - 1,
                    });
                }
                if ui.button("再読込").clicked() {
                    let _ = self.sender.send(WorkerCommand::SelectFocused {
                        effect_index: effect_number - 1,
                    });
                }
                ui.separator();
                ui.label("追加");
                egui::ComboBox::from_id_salt("puppet_pin_kind")
                    .selected_text(pin_kind_text(self.pin_kind, ui.visuals().text_color()))
                    .show_ui(ui, |ui| {
                        for kind in PinKind::ALL {
                            let text = pin_kind_text(kind, ui.visuals().text_color());
                            ui.selectable_value(&mut self.pin_kind, kind, text);
                        }
                    });
                ui.separator();
                if ui.button("全体表示").clicked() {
                    self.view = ViewState::default();
                }
            });
        });
    }

    fn canvas(&mut self, ui: &mut egui::Ui, frame: &CapturedFrame) {
        let editor = self.shared.editor_snapshot();
        let available = ui.available_size().max(egui::vec2(1.0, 1.0));
        let (response, painter) = ui.allocate_painter(available, egui::Sense::click_and_drag());
        let canvas = response.rect;
        painter.rect_filled(canvas, 0.0, ui.visuals().extreme_bg_color);

        let image_size = egui::vec2(frame.width as f32, frame.height as f32);
        let fit = ((canvas.width() - IMAGE_MARGIN * 2.0) / image_size.x)
            .min((canvas.height() - IMAGE_MARGIN * 2.0) / image_size.y)
            .max(0.001);
        let scale = fit * self.view.zoom;
        let displayed_size = image_size * scale;
        let image_rect =
            egui::Rect::from_center_size(canvas.center() + self.view.pan, displayed_size);
        self.help_tooltip(ui, &response, image_rect);

        draw_checkerboard(&painter, image_rect.intersect(canvas), scale);
        if let Some(texture) = &self.texture {
            painter.image(
                texture.id(),
                image_rect,
                egui::Rect::from_min_max(egui::Pos2::ZERO, egui::pos2(1.0, 1.0)),
                egui::Color32::WHITE,
            );
        }

        let pointer = response.hover_pos();
        let hovered_pin = pointer.and_then(|position| {
            self.local_model
                .as_ref()
                .and_then(|model| nearest_pin(model, position, image_rect, scale, PIN_HOVER_RADIUS))
        });
        let hovered_bone_connection = pointer.and_then(|position| {
            (hovered_pin.is_none()
                && self
                    .local_model
                    .as_ref()
                    .is_some_and(|model| model.src.len() < MAX_PINS))
            .then(|| {
                self.local_model.as_ref().and_then(|model| {
                    nearest_bone_connection(
                        model,
                        position,
                        image_rect,
                        scale,
                        BONE_LINE_HOVER_RADIUS,
                    )
                })
            })
            .flatten()
        });
        if let Some(model) = &self.local_model {
            for (child, parent) in model.bone_parents.iter().copied().enumerate() {
                if let Some(parent) = parent {
                    let from = point_to_screen(model.src[parent], image_rect, scale);
                    let to = point_to_screen(model.src[child], image_rect, scale);
                    let hovered = hovered_bone_connection
                        .is_some_and(|hit| hit.parent == parent && hit.child == child);
                    draw_bone_connection(&painter, from, to, hovered);
                }
            }
            if let Some(hit) = hovered_bone_connection {
                draw_pin_shape(
                    &painter,
                    hit.position,
                    PIN_RADIUS * 1.35,
                    PinKind::Bone,
                    pin_color(PinKind::Bone).gamma_multiply(0.38),
                );
            }
            for (index, point) in model.src.iter().enumerate() {
                let position = point_to_screen(*point, image_rect, scale);
                let hovered = hovered_pin == Some(index) || self.dragging_pin == Some(index);
                let radius = if hovered {
                    PIN_RADIUS * 1.65
                } else {
                    PIN_RADIUS
                };
                draw_pin_shape(
                    &painter,
                    position,
                    radius,
                    model.kinds[index],
                    pin_color(model.kinds[index]),
                );
                if self.selected_pin == Some(index) {
                    painter.circle_stroke(
                        position,
                        PIN_RADIUS + 3.0,
                        egui::Stroke::new(1.5, pin_color(model.kinds[index])),
                    );
                }
                painter.text(
                    position + egui::vec2(-radius - 2.0, radius + 1.0),
                    egui::Align2::RIGHT_TOP,
                    format!("#{}", index + 1),
                    egui::FontId::monospace(10.0),
                    egui::Color32::WHITE,
                );
            }
        }

        if response.hovered() {
            let (scroll, zoom_delta, modifiers, pointer_position, middle_down, pointer_delta) = ui
                .input(|input| {
                    (
                        input.smooth_scroll_delta(),
                        input.zoom_delta(),
                        input.modifiers,
                        input.pointer.hover_pos(),
                        input.pointer.button_down(egui::PointerButton::Middle),
                        input.pointer.delta(),
                    )
                });
            if modifiers.ctrl && zoom_delta != 1.0 {
                zoom_about_pointer(
                    &mut self.view,
                    zoom_delta,
                    pointer_position,
                    canvas.center(),
                );
            } else {
                self.view.pan += wheel_pan_delta(scroll, modifiers.shift);
            }
            if middle_down {
                self.view.pan += pointer_delta;
            }
        }

        if response.drag_started_by(egui::PointerButton::Primary)
            && !editor.busy
            && let Some(index) = hovered_pin
        {
            self.dragging_pin = Some(index);
        }

        let primary_double_clicked = response.double_clicked_by(egui::PointerButton::Primary);
        if primary_double_clicked && !editor.busy {
            if let Some(index) = hovered_pin {
                self.selected_pin = Some(index);
                if let Some(kind) = self
                    .local_model
                    .as_ref()
                    .and_then(|model| model.kinds.get(index))
                    .copied()
                    && kind != PinKind::Starch
                    && let Some(target) = editor.target.clone()
                {
                    let _ = self.sender.send(WorkerCommand::AddPinControl {
                        target,
                        pin_index: index,
                        kind,
                    });
                }
            }
        } else if response.clicked_by(egui::PointerButton::Primary) && !editor.busy {
            if let Some(index) = hovered_pin {
                let shift = ui.input(|input| input.modifiers.shift);
                let attachment = shift.then(|| {
                    let model = self.local_model.as_ref()?;
                    let parent = self.selected_pin?;
                    (parent != index
                        && model.kinds.get(parent) == Some(&PinKind::Bone)
                        && model.kinds.get(index) == Some(&PinKind::Bone))
                    .then_some(PinOperation::AttachRoot {
                        parent,
                        root: index,
                    })
                });
                if let Some(operation) = attachment.flatten() {
                    self.optimistic_apply(operation);
                } else {
                    self.selected_pin = Some(index);
                }
            } else if let Some(hit) = hovered_bone_connection {
                let point = screen_to_point(hit.position, image_rect, scale);
                self.optimistic_apply(PinOperation::InsertBone {
                    point,
                    parent: hit.parent,
                    child: hit.child,
                });
            } else if let Some(position) = response.interact_pointer_pos()
                && image_rect.contains(position)
                && self
                    .local_model
                    .as_ref()
                    .is_some_and(|model| model.src.len() < MAX_PINS)
            {
                let bone_parent = (self.pin_kind == PinKind::Bone)
                    .then(|| {
                        self.selected_pin.filter(|&selected| {
                            self.local_model
                                .as_ref()
                                .and_then(|model| model.kinds.get(selected))
                                == Some(&PinKind::Bone)
                        })
                    })
                    .flatten();
                let point = screen_to_point(position, image_rect, scale);
                self.optimistic_apply(PinOperation::Add {
                    point,
                    kind: self.pin_kind,
                    bone_parent,
                });
            }
        }

        if let Some(index) = self.dragging_pin {
            let (primary_down, released, pointer_position) = ui.input(|input| {
                (
                    input.pointer.button_down(egui::PointerButton::Primary),
                    input.pointer.button_released(egui::PointerButton::Primary),
                    input.pointer.interact_pos(),
                )
            });
            if primary_down
                && let Some(position) = pointer_position
                && let Some(model) = &mut self.local_model
                && let Some(point) = model.src.get_mut(index)
            {
                *point = screen_to_point(position, image_rect, scale);
            }
            if released || !primary_down {
                if let Some(point) = self
                    .local_model
                    .as_ref()
                    .and_then(|model| model.src.get(index))
                    .copied()
                {
                    self.send_operation(PinOperation::MoveSourceAndDestination { index, point });
                }
                self.dragging_pin = None;
            }
        }

        if response.clicked_by(egui::PointerButton::Secondary) && !editor.busy {
            if let Some(index) = hovered_pin
                && self
                    .local_model
                    .as_ref()
                    .is_some_and(|model| model.src.len() > MIN_PINS)
            {
                self.optimistic_apply(PinOperation::Delete(index));
            } else if hovered_pin.is_none() {
                self.selected_pin = None;
            }
        }
    }

    fn help_tooltip(
        &mut self,
        ui: &mut egui::Ui,
        response: &egui::Response,
        image_rect: egui::Rect,
    ) {
        let pointer = response.hover_pos();
        let now = ui.input(|input| input.time);
        let moved = pointer != self.last_help_pointer;
        self.last_help_pointer = pointer;

        let eligible = pointer.is_some_and(|position| !image_rect.contains(position))
            && !ui.input(|input| input.pointer.any_down());
        if !eligible {
            self.help_hover_started_at = None;
            return;
        }
        if moved || self.help_hover_started_at.is_none() {
            self.help_hover_started_at = Some(now);
        }
        let elapsed = now - self.help_hover_started_at.unwrap_or(now);
        if elapsed >= HELP_TOOLTIP_DELAY_SECONDS {
            egui::Tooltip::for_widget(response).at_pointer().show(|ui| {
                ui.label(HELP_TEXT);
            });
        } else {
            ui.ctx()
                .request_repaint_after(std::time::Duration::from_secs_f64(
                    HELP_TOOLTIP_DELAY_SECONDS - elapsed,
                ));
        }
    }

    fn optimistic_apply(&mut self, operation: PinOperation) {
        if let Some(model) = &mut self.local_model
            && model.apply(operation).is_ok()
        {
            match operation {
                PinOperation::Add { .. } | PinOperation::InsertBone { .. } => {
                    self.selected_pin = Some(model.src.len() - 1);
                }
                PinOperation::Delete(deleted) => {
                    self.selected_pin = match self.selected_pin {
                        Some(selected) if selected == deleted => None,
                        Some(selected) if selected > deleted => Some(selected - 1),
                        selected => selected,
                    };
                }
                PinOperation::AttachRoot { root, .. } => self.selected_pin = Some(root),
                PinOperation::MoveSourceAndDestination { .. } => {}
            }
            self.send_operation(operation);
        }
    }

    fn send_operation(&self, operation: PinOperation) {
        let editor = self.shared.editor_snapshot();
        if let Some(target) = editor.target {
            let _ = self.sender.send(WorkerCommand::Apply { target, operation });
        }
    }
}

fn pin_kind_text(kind: PinKind, text_color: egui::Color32) -> egui::text::LayoutJob {
    let symbol = match kind {
        PinKind::Bone => "◆ ",
        PinKind::Starch | PinKind::Overlap => "○ ",
        _ => "● ",
    };
    let mut job = egui::text::LayoutJob::default();
    job.append(
        symbol,
        0.0,
        egui::TextFormat {
            color: pin_color(kind),
            ..Default::default()
        },
    );
    job.append(
        kind.label(),
        0.0,
        egui::TextFormat {
            color: text_color,
            ..Default::default()
        },
    );
    job
}

fn pin_color(kind: PinKind) -> egui::Color32 {
    match kind {
        PinKind::Position => egui::Color32::from_rgb(242, 191, 38),
        PinKind::Bone => egui::Color32::from_rgb(32, 224, 96),
        PinKind::Bend => egui::Color32::from_rgb(255, 160, 32),
        PinKind::Detail => egui::Color32::from_rgb(224, 40, 72),
        PinKind::Starch => egui::Color32::from_rgb(160, 96, 224),
        PinKind::Overlap => egui::Color32::from_rgb(64, 128, 255),
    }
}

fn draw_pin_shape(
    painter: &egui::Painter,
    position: egui::Pos2,
    radius: f32,
    kind: PinKind,
    color: egui::Color32,
) {
    let outline = if color.a() == 255 {
        egui::Color32::BLACK
    } else {
        color
    };
    match kind {
        PinKind::Bone => {
            painter.add(egui::Shape::convex_polygon(
                vec![
                    position + egui::vec2(0.0, -radius),
                    position + egui::vec2(radius, 0.0),
                    position + egui::vec2(0.0, radius),
                    position + egui::vec2(-radius, 0.0),
                ],
                color,
                egui::Stroke::new(1.5, outline),
            ));
        }
        PinKind::Starch | PinKind::Overlap => {
            let outline_width = 4.5;
            let color_width = 3.5;
            let stroke_radius = (radius - outline_width * 0.5).max(0.5);
            if color.a() == 255 {
                painter.circle_stroke(
                    position,
                    stroke_radius,
                    egui::Stroke::new(outline_width, outline),
                );
            }
            painter.circle_stroke(
                position,
                stroke_radius,
                egui::Stroke::new(color_width, color),
            );
        }
        _ => {
            painter.circle_filled(position, radius, color);
            painter.circle_stroke(position, radius, egui::Stroke::new(1.5, outline));
        }
    }
}

fn draw_bone_connection(painter: &egui::Painter, from: egui::Pos2, to: egui::Pos2, hovered: bool) {
    let delta = to - from;
    let length = delta.length();
    if length <= f32::EPSILON {
        return;
    }
    let direction = delta / length;
    let normal = egui::vec2(-direction.y, direction.x);
    let color = if hovered {
        pin_color(PinKind::Bone)
    } else {
        egui::Color32::from_rgb(96, 220, 128)
    };
    painter.line_segment(
        [from, to],
        egui::Stroke::new(if hovered { 3.0 } else { 2.0 }, color),
    );

    let center = from + delta * 0.58;
    let tip = center + direction * BONE_ARROW_SIZE;
    let base = center - direction * BONE_ARROW_SIZE * 0.75;
    painter.add(egui::Shape::convex_polygon(
        vec![
            tip,
            base + normal * BONE_ARROW_SIZE * 0.65,
            base - normal * BONE_ARROW_SIZE * 0.65,
        ],
        color,
        egui::Stroke::NONE,
    ));
}

impl eframe::App for PuppetEditorApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.refresh_local_state();
        self.toolbar(ui);
        let editor = self.shared.editor_snapshot();
        let capture_status = self.shared.capture_snapshot().status;
        egui::Panel::bottom("puppet_status").show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(&editor.status);
                ui.separator();
                ui.weak(capture_status);
            });
        });

        let frame = self.refresh_texture(ui.ctx());
        egui::CentralPanel::default().show(ui, |ui| {
            if editor.target.is_none() {
                ui.centered_and_justified(|ui| ui.label(&editor.status));
            } else if editor.model.is_none() {
                ui.centered_and_justified(|ui| ui.label(&editor.status));
            } else if let Some(frame) = frame {
                self.canvas(ui, &frame);
            } else {
                ui.centered_and_justified(|ui| ui.label("変形前画像を待機しています。"));
            }
        });
    }
}

fn point_to_screen(point: Point, image_rect: egui::Rect, scale: f32) -> egui::Pos2 {
    egui::pos2(
        image_rect.center().x + point.x * scale,
        image_rect.center().y + point.y * scale,
    )
}

fn screen_to_point(position: egui::Pos2, image_rect: egui::Rect, scale: f32) -> Point {
    Point {
        x: (position.x - image_rect.center().x) / scale,
        y: (position.y - image_rect.center().y) / scale,
    }
}

fn nearest_pin(
    model: &PinModel,
    pointer: egui::Pos2,
    image_rect: egui::Rect,
    scale: f32,
    radius: f32,
) -> Option<usize> {
    model
        .src
        .iter()
        .enumerate()
        .filter_map(|(index, point)| {
            let distance = point_to_screen(*point, image_rect, scale).distance(pointer);
            (distance <= radius).then_some((index, distance))
        })
        .min_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(index, _)| index)
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct BoneConnectionHit {
    parent: usize,
    child: usize,
    position: egui::Pos2,
}

fn nearest_bone_connection(
    model: &PinModel,
    pointer: egui::Pos2,
    image_rect: egui::Rect,
    scale: f32,
    radius: f32,
) -> Option<BoneConnectionHit> {
    model
        .bone_parents
        .iter()
        .copied()
        .enumerate()
        .filter_map(|(child, parent)| {
            let parent = parent?;
            let from = point_to_screen(model.src[parent], image_rect, scale);
            let to = point_to_screen(model.src[child], image_rect, scale);
            let segment = to - from;
            let length_squared = segment.length_sq();
            if length_squared <= f32::EPSILON {
                return None;
            }
            let offset = pointer - from;
            let projection = (offset.dot(segment) / length_squared).clamp(0.0, 1.0);
            let position = from + segment * projection;
            let distance = position.distance(pointer);
            (distance <= radius).then_some((
                BoneConnectionHit {
                    parent,
                    child,
                    position,
                },
                distance,
            ))
        })
        .min_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(hit, _)| hit)
}

fn wheel_pan_delta(scroll: egui::Vec2, shift: bool) -> egui::Vec2 {
    if shift {
        egui::vec2(scroll.x + scroll.y, 0.0)
    } else {
        egui::vec2(0.0, scroll.y)
    }
}

fn zoom_about_pointer(
    view: &mut ViewState,
    zoom_delta: f32,
    pointer: Option<egui::Pos2>,
    canvas_center: egui::Pos2,
) {
    let old_zoom = view.zoom;
    view.zoom = (view.zoom * zoom_delta).clamp(MIN_ZOOM, MAX_ZOOM);
    if let Some(pointer) = pointer {
        let image_center = canvas_center + view.pan;
        let ratio = view.zoom / old_zoom;
        view.pan += (pointer - image_center) * (1.0 - ratio);
    }
}

fn draw_checkerboard(painter: &egui::Painter, rect: egui::Rect, scale: f32) {
    let tile = (12.0 * scale.sqrt()).clamp(6.0, 20.0);
    let light = egui::Color32::from_gray(92);
    let dark = egui::Color32::from_gray(62);
    let cols = (rect.width() / tile).ceil() as usize;
    let rows = (rect.height() / tile).ceil() as usize;
    for y in 0..rows {
        for x in 0..cols {
            let min = rect.min + egui::vec2(x as f32 * tile, y as f32 * tile);
            let tile_rect = egui::Rect::from_min_size(min, egui::vec2(tile, tile)).intersect(rect);
            painter.rect_filled(tile_rect, 0.0, if (x + y) % 2 == 0 { light } else { dark });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn image_and_script_coordinates_round_trip() {
        let rect = egui::Rect::from_center_size(egui::pos2(200.0, 150.0), egui::vec2(400.0, 200.0));
        let point = Point { x: -12.5, y: 30.0 };
        let screen = point_to_screen(point, rect, 2.0);
        assert_eq!(screen_to_point(screen, rect, 2.0), point);
    }

    #[test]
    fn pin_hit_radius_is_in_screen_pixels() {
        let model = PinModel::parse("1", "0,0", "0,0", "", "").unwrap();
        let rect = egui::Rect::from_center_size(egui::pos2(100.0, 100.0), egui::vec2(100.0, 100.0));
        assert_eq!(
            nearest_pin(&model, egui::pos2(111.0, 100.0), rect, 0.2, 12.0),
            Some(0)
        );
        assert_eq!(
            nearest_pin(&model, egui::pos2(113.0, 100.0), rect, 10.0, 12.0),
            None
        );
    }

    #[test]
    fn bone_connection_hit_returns_the_projected_insertion_point() {
        let model = PinModel::parse("2", "-50,0,50,0", "-50,0,50,0", "1,1", "{{1,{2}}}").unwrap();
        let rect = egui::Rect::from_center_size(egui::pos2(100.0, 100.0), egui::vec2(200.0, 200.0));
        let hit =
            nearest_bone_connection(&model, egui::pos2(125.0, 106.0), rect, 1.0, 10.0).unwrap();
        assert_eq!((hit.parent, hit.child), (0, 1));
        assert_eq!(hit.position, egui::pos2(125.0, 100.0));
        assert!(
            nearest_bone_connection(&model, egui::pos2(125.0, 111.0), rect, 1.0, 10.0).is_none()
        );
    }

    #[test]
    fn shift_wheel_moves_horizontally_for_transformed_or_raw_delta() {
        assert_eq!(
            wheel_pan_delta(egui::vec2(40.0, 0.0), true),
            egui::vec2(40.0, 0.0)
        );
        assert_eq!(
            wheel_pan_delta(egui::vec2(0.0, 40.0), true),
            egui::vec2(40.0, 0.0)
        );
        assert_eq!(
            wheel_pan_delta(egui::vec2(40.0, 20.0), false),
            egui::vec2(0.0, 20.0)
        );
    }

    #[test]
    fn zoom_keeps_the_pointer_over_the_same_image_position() {
        let mut view = ViewState::default();
        let center = egui::pos2(100.0, 100.0);
        let pointer = egui::pos2(150.0, 125.0);
        zoom_about_pointer(&mut view, 2.0, Some(pointer), center);
        assert_eq!(view.zoom, 2.0);
        assert_eq!(view.pan, egui::vec2(-50.0, -25.0));
    }
}
