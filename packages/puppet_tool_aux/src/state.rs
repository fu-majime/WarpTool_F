use std::sync::{
    Arc, Mutex, OnceLock,
    atomic::{AtomicU64, Ordering},
};

use aviutl2_eframe::egui;

use crate::model::PinModel;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EditTarget {
    pub object: aviutl2::generic::ObjectHandle,
    pub effect_name: String,
    pub effect_index: usize,
}

#[derive(Clone, Debug)]
pub struct EditorSnapshot {
    pub target: Option<EditTarget>,
    pub model: Option<PinModel>,
    pub status: String,
    pub busy: bool,
    pub revision: u64,
    pub capture_after: u64,
}

impl Default for EditorSnapshot {
    fn default() -> Self {
        Self {
            target: None,
            model: None,
            status: "タイムラインで編集するオブジェクトを選択してください。".to_owned(),
            busy: false,
            revision: 0,
            capture_after: 0,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CaptureSource {
    InputImage,
    FramebufferFallback,
}

#[derive(Clone, Debug)]
pub struct CapturedFrame {
    pub sequence: u64,
    pub width: usize,
    pub height: usize,
    pub rgba: Arc<[u8]>,
}

#[derive(Clone, Debug)]
pub struct CaptureSnapshot {
    pub latest: Option<CapturedFrame>,
    pub status: String,
}

impl Default for CaptureSnapshot {
    fn default() -> Self {
        Self {
            latest: None,
            status: "変形前画像を待機しています。".to_owned(),
        }
    }
}

pub struct SharedState {
    editor: Mutex<EditorSnapshot>,
    capture: Mutex<CaptureSnapshot>,
    capture_sequence: AtomicU64,
    egui_context: OnceLock<egui::Context>,
}

impl SharedState {
    pub fn new() -> Self {
        Self {
            editor: Mutex::new(EditorSnapshot::default()),
            capture: Mutex::new(CaptureSnapshot::default()),
            capture_sequence: AtomicU64::new(0),
            egui_context: OnceLock::new(),
        }
    }

    pub fn set_egui_context(&self, context: egui::Context) {
        let _ = self.egui_context.set(context);
    }

    pub fn repaint(&self) {
        if let Some(context) = self.egui_context.get() {
            context.request_repaint();
        }
    }

    pub fn editor_snapshot(&self) -> EditorSnapshot {
        self.editor.lock().unwrap().clone()
    }

    pub fn update_editor(&self, update: impl FnOnce(&mut EditorSnapshot)) {
        let mut editor = self.editor.lock().unwrap();
        update(&mut editor);
        editor.revision = editor.revision.wrapping_add(1);
        drop(editor);
        self.repaint();
    }

    pub fn capture_snapshot(&self) -> CaptureSnapshot {
        self.capture.lock().unwrap().clone()
    }

    pub fn current_capture_sequence(&self) -> u64 {
        self.capture_sequence.load(Ordering::Acquire)
    }

    pub fn publish_capture(
        &self,
        width: usize,
        height: usize,
        rgba: Vec<u8>,
        source: CaptureSource,
    ) {
        let sequence = self.capture_sequence.fetch_add(1, Ordering::AcqRel) + 1;
        let source_label = match source {
            CaptureSource::InputImage => "入力画像",
            CaptureSource::FramebufferFallback => "フレームバッファ（フォールバック）",
        };
        let mut capture = self.capture.lock().unwrap();
        capture.status = format!("{source_label}: {width} x {height}");
        capture.latest = Some(CapturedFrame {
            sequence,
            width,
            height,
            rgba: rgba.into(),
        });
        drop(capture);
        self.repaint();
    }

    pub fn publish_capture_error(&self, status: impl Into<String>) {
        let mut capture = self.capture.lock().unwrap();
        capture.status = status.into();
        drop(capture);
        self.repaint();
    }
}

static SHARED_STATE: OnceLock<Arc<SharedState>> = OnceLock::new();

pub fn shared_state() -> Arc<SharedState> {
    SHARED_STATE
        .get_or_init(|| Arc::new(SharedState::new()))
        .clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capture_publication_advances_the_generation() {
        let shared = SharedState::new();
        assert_eq!(shared.current_capture_sequence(), 0);
        shared.publish_capture(1, 1, vec![1, 2, 3, 4], CaptureSource::InputImage);
        let snapshot = shared.capture_snapshot();
        let first = snapshot.latest.unwrap();
        assert_eq!(first.sequence, 1);
        assert!(snapshot.status.contains("入力画像"));
        shared.publish_capture(1, 1, vec![5, 6, 7, 8], CaptureSource::FramebufferFallback);
        assert_eq!(shared.capture_snapshot().latest.unwrap().sequence, 2);
    }
}
