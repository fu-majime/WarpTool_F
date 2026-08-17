//! Minimal bridge from an AviUtl2 Lua script to the project editing API.
//!
//! This DLL is loaded as a generic plugin (`.aux2`). It registers the
//! `ScriptEditBridge` script module internally. Script calls only enqueue an
//! update; a worker thread performs `call_edit_section` and
//! `set_object_item_value`, avoiding project edits from inside Lua rendering.

use std::{
    sync::{
        OnceLock,
        mpsc::{Receiver, SyncSender, TrySendError, sync_channel},
    },
    thread::JoinHandle,
};

use aviutl2::{AnyResult, module::ScriptModuleFunctions};

const MODULE_NAME: &str = "ScriptEditBridge";
const MAX_EFFECT_INDEX: i32 = 63;

static EDIT_HANDLE: aviutl2::generic::GlobalEditHandle = aviutl2::generic::GlobalEditHandle::new();
static REQUEST_SENDER: OnceLock<SyncSender<WorkerMessage>> = OnceLock::new();

#[derive(Debug)]
struct ItemRequest {
    effect_name: String,
    effect_index: usize,
    item_name: String,
    value: String,
}

#[derive(Debug)]
enum WorkerMessage {
    SetItem(ItemRequest),
    Shutdown,
}

#[aviutl2::plugin(ScriptModule)]
struct ScriptEditModule;

impl aviutl2::module::ScriptModule for ScriptEditModule {
    fn new(_info: aviutl2::AviUtl2Info) -> AnyResult<Self> {
        Ok(Self)
    }

    fn plugin_info(&self) -> aviutl2::module::ScriptModuleTable {
        aviutl2::module::ScriptModuleTable {
            information: format!(
                "Script Lua/edit bridge / v{} / Fu-Majime",
                env!("CARGO_PKG_VERSION")
            ),
            functions: Self::functions(),
        }
    }
}

#[aviutl2::module::functions]
impl ScriptEditModule {
    /// Queue a replacement for an effect item's value on the focused object.
    ///
    /// `false` means the bridge is unavailable, an argument is invalid, or an
    /// update is already waiting. Lua rendering must never wait for the
    /// project-editing thread.
    fn request(effect_name: String, effect_index: i32, item_name: String, value: String) -> bool {
        if effect_name.is_empty()
            || item_name.is_empty()
            || !(0..=MAX_EFFECT_INDEX).contains(&effect_index)
        {
            return false;
        }

        let Some(sender) = REQUEST_SENDER.get() else {
            return false;
        };
        let request = ItemRequest {
            effect_name,
            effect_index: effect_index as usize,
            item_name,
            value,
        };

        match sender.try_send(WorkerMessage::SetItem(request)) {
            Ok(()) => true,
            Err(TrySendError::Full(_) | TrySendError::Disconnected(_)) => false,
        }
    }
}

#[aviutl2::plugin(GenericPlugin)]
struct ScriptEditBridgePlugin {
    script_module: aviutl2::generic::SubPlugin<ScriptEditModule>,
    worker: Option<JoinHandle<()>>,
}

impl aviutl2::generic::GenericPlugin for ScriptEditBridgePlugin {
    fn new(info: aviutl2::AviUtl2Info) -> AnyResult<Self> {
        let script_module = aviutl2::generic::SubPlugin::new_script_module(&info)?;

        // A capacity of one coalesces the many calls a render script can make.
        // A later render retries if the queue was full.
        let (sender, receiver) = sync_channel(1);
        REQUEST_SENDER
            .set(sender)
            .map_err(|_| aviutl2::anyhow::anyhow!("request queue was already initialized"))?;
        let worker = std::thread::Builder::new()
            .name("script-edit-worker".to_owned())
            .spawn(move || worker_loop(receiver))?;

        Ok(Self {
            script_module,
            worker: Some(worker),
        })
    }

    fn register(&mut self, registry: &mut aviutl2::generic::HostAppHandle) {
        EDIT_HANDLE.init(registry.create_edit_handle());
        registry.register_script_module(Some(MODULE_NAME), &self.script_module);
    }

    fn plugin_info(&self) -> aviutl2::generic::GenericPluginTable {
        aviutl2::generic::GenericPluginTable {
            name: "Script Edit Bridge".to_owned(),
            information: format!(
                "Queued Lua-to-edit bridge / v{} / Fu-Majime",
                env!("CARGO_PKG_VERSION")
            ),
        }
    }
}

impl Drop for ScriptEditBridgePlugin {
    fn drop(&mut self) {
        if let Some(sender) = REQUEST_SENDER.get() {
            // Blocking here is intentional: if one request is pending, let the
            // worker finish it before asking the worker to stop.
            let _ = sender.send(WorkerMessage::Shutdown);
        }
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

fn worker_loop(receiver: Receiver<WorkerMessage>) {
    while let Ok(message) = receiver.recv() {
        match message {
            WorkerMessage::SetItem(request) => apply_item(request),
            WorkerMessage::Shutdown => break,
        }
    }
}

fn apply_item(request: ItemRequest) {
    if !EDIT_HANDLE.is_ready() {
        return;
    }

    let result = EDIT_HANDLE.call_edit_section(move |edit_section| -> AnyResult<bool> {
        let Some(object) = edit_section.get_focused_object()? else {
            return Ok(false);
        };

        // Reading first both verifies that the focused object owns the item
        // and avoids creating needless undo entries for an unchanged value.
        let Ok(current) = edit_section.get_object_effect_item(
            object,
            &request.effect_name,
            request.effect_index,
            &request.item_name,
        ) else {
            return Ok(false);
        };
        if current == request.value {
            return Ok(false);
        }

        edit_section.set_object_effect_item(
            object,
            &request.effect_name,
            request.effect_index,
            &request.item_name,
            &request.value,
        )?;
        Ok(true)
    });

    match result {
        Ok(Ok(_)) => {}
        Ok(Err(error)) => aviutl2::tracing::warn!("Script edit update failed: {error:#}"),
        Err(error) => aviutl2::tracing::warn!("call_edit_section failed: {error}"),
    }
}

aviutl2::register_generic_plugin!(ScriptEditBridgePlugin);
