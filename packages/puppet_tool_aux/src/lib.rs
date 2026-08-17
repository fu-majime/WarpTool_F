//! Native editor and input-frame capture bridge for WarpTool_F's puppet script.

use std::{
    sync::{
        Arc,
        mpsc::{Sender, channel},
    },
    thread::JoinHandle,
};

use aviutl2::AnyResult;

mod capture;
mod model;
mod state;
mod ui;
mod worker;

use capture::FrameCapture;
use state::{SharedState, shared_state};
use worker::WorkerCommand;

pub const INTERNAL_CAPTURE_FILTER_NAME: &str = "WarpTool_F Internal Puppet Capture";
pub const PUPPET_EFFECT_NAME: &str = "パペットツール@WarpTool_F";
const WINDOW_NAME: &str = "WarpTool_F パペットツール";

pub static EDIT_HANDLE: aviutl2::generic::GlobalEditHandle =
    aviutl2::generic::GlobalEditHandle::new();

#[aviutl2::plugin(FilterPlugin)]
struct PuppetCaptureFilter {
    capture: FrameCapture,
    shared: Arc<SharedState>,
}

impl aviutl2::filter::FilterPlugin for PuppetCaptureFilter {
    type Userdata = ();

    fn new(_info: aviutl2::AviUtl2Info) -> AnyResult<Self> {
        Ok(Self {
            capture: FrameCapture::new(),
            shared: shared_state(),
        })
    }

    fn plugin_info(&self) -> aviutl2::filter::FilterPluginTable {
        aviutl2::filter::FilterPluginTable {
            name: INTERNAL_CAPTURE_FILTER_NAME.to_owned(),
            label: None,
            information: format!(
                "WarpTool_F internal input capture / v{} / Fu-Majime",
                env!("CARGO_PKG_VERSION")
            ),
            flags: aviutl2::bitflag!(aviutl2::filter::FilterPluginFlags { video: true }),
            config_items: Vec::new(),
        }
    }

    fn proc_video(
        &self,
        _config: &[aviutl2::filter::FilterConfigItem],
        video: &mut aviutl2::filter::FilterProcVideo<Self::Userdata>,
    ) -> AnyResult<()> {
        self.capture.capture(video, &self.shared);
        Ok(())
    }
}

#[aviutl2::plugin(GenericPlugin)]
struct PuppetToolAux {
    window: aviutl2_eframe::EframeWindow,
    capture_filter: aviutl2::generic::SubPlugin<PuppetCaptureFilter>,
    sender: Sender<WorkerCommand>,
    worker: Option<JoinHandle<()>>,
}

impl aviutl2::generic::GenericPlugin for PuppetToolAux {
    fn new(info: aviutl2::AviUtl2Info) -> AnyResult<Self> {
        let shared = shared_state();
        let (sender, receiver) = channel();
        let worker_shared = shared.clone();
        let worker = std::thread::Builder::new()
            .name("warp-tool-puppet-editor".to_owned())
            .spawn(move || worker::run(receiver, worker_shared))?;

        let app_shared = shared.clone();
        let app_sender = sender.clone();
        let window = aviutl2_eframe::EframeWindow::new(
            "WarpToolFPuppetEditor",
            move |creation_context, handle| {
                Ok(Box::new(ui::PuppetEditorApp::new(
                    creation_context,
                    handle,
                    app_shared.clone(),
                    app_sender.clone(),
                )))
            },
        )?;

        Ok(Self {
            window,
            capture_filter: aviutl2::generic::SubPlugin::new_filter_plugin(&info)?,
            sender,
            worker: Some(worker),
        })
    }

    fn plugin_info(&self) -> aviutl2::generic::GenericPluginTable {
        aviutl2::generic::GenericPluginTable {
            name: "WarpTool_F Puppet Editor".to_owned(),
            information: format!(
                "WarpTool_F puppet pin editor / v{} / Fu-Majime",
                env!("CARGO_PKG_VERSION")
            ),
        }
    }

    fn register(&mut self, registry: &mut aviutl2::generic::HostAppHandle) {
        EDIT_HANDLE.init(registry.create_edit_handle());
        registry.register_filter_plugin(&self.capture_filter);
        if let Ok(handle) = self.window.handle()
            && let Err(error) = registry.register_window_client(WINDOW_NAME, &handle)
        {
            aviutl2::tracing::error!("failed to register puppet editor window: {error}");
        }
    }

    fn event_update_object_info(&mut self) {
        let _ = self.sender.send(WorkerCommand::Refresh);
    }

    fn event_change_focus_object(&mut self) {
        let _ = self
            .sender
            .send(WorkerCommand::SelectFocused { effect_index: 0 });
    }

    fn on_project_load(&mut self, _project: &mut aviutl2::generic::ProjectFile) {
        let _ = self
            .sender
            .send(WorkerCommand::SelectFocused { effect_index: 0 });
    }
}

impl Drop for PuppetToolAux {
    fn drop(&mut self) {
        let _ = self.sender.send(WorkerCommand::Shutdown);
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

aviutl2::register_generic_plugin!(PuppetToolAux);
