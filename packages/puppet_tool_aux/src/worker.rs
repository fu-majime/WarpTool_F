use std::sync::{Arc, mpsc::Receiver};

use aviutl2::AnyResult;
use aviutl2::anyhow::Context as _;
use aviutl2::generic::ReadSectionProvider;

use crate::{
    EDIT_HANDLE, PUPPET_EFFECT_NAME,
    model::{PinKind, PinModel, PinOperation, format_number},
    state::{EditTarget, SharedState},
};

// AviUtl2's edit API addresses script settings by their alias/display names,
// not by the Lua variable names after `@`.
const PINS_ITEM_NAME: &str = "ピン数";
const SOURCE_ITEM_NAME: &str = "ピン元";
const DESTINATION_ITEM_NAME: &str = "ピン先";
const PIN_TYPES_ITEM_NAME: &str = "ピン種類";
const BONE_FOREST_ITEM_NAME: &str = "ボーン親子関係";
const CONTROL_PIN_NUMBER_ITEM_NAME: &str = "番号";
const CONTROL_X_ITEM_NAME: &str = "X";
const CONTROL_Y_ITEM_NAME: &str = "Y";
const POSITION_CONTROL_EFFECT_NAME: &str = "パペットツール@位置ピン制御@WarpTool_F";
const BONE_CONTROL_EFFECT_NAME: &str = "パペットツール@ボーンピン制御@WarpTool_F";
const BEND_CONTROL_EFFECT_NAME: &str = "パペットツール@ベンドピン制御@WarpTool_F";
const DETAIL_CONTROL_EFFECT_NAME: &str = "パペットツール@詳細ピン制御@WarpTool_F";
const OVERLAP_CONTROL_EFFECT_NAME: &str = "パペットツール@重なりピン制御@WarpTool_F";
const NEW_OBJECT_LENGTH: usize = 81;

#[derive(Clone, Debug)]
pub enum WorkerCommand {
    SelectFocused {
        effect_index: usize,
    },
    Refresh,
    Apply {
        target: EditTarget,
        operation: PinOperation,
    },
    AddPinControl {
        target: EditTarget,
        pin_index: usize,
        kind: PinKind,
    },
    AddObject,
    Shutdown,
}

pub fn run(receiver: Receiver<WorkerCommand>, shared: Arc<SharedState>) {
    while let Ok(command) = receiver.recv() {
        match command {
            WorkerCommand::SelectFocused { effect_index } => select_focused(&shared, effect_index),
            WorkerCommand::Refresh => refresh(&shared),
            WorkerCommand::Apply { target, operation } => apply(&shared, target, operation),
            WorkerCommand::AddPinControl {
                target,
                pin_index,
                kind,
            } => add_pin_control(&shared, target, pin_index, kind),
            WorkerCommand::AddObject => add_object(&shared),
            WorkerCommand::Shutdown => break,
        }
    }
}

fn control_effect_name(kind: PinKind) -> Option<&'static str> {
    match kind {
        PinKind::Position => Some(POSITION_CONTROL_EFFECT_NAME),
        PinKind::Bone => Some(BONE_CONTROL_EFFECT_NAME),
        PinKind::Bend => Some(BEND_CONTROL_EFFECT_NAME),
        PinKind::Detail => Some(DETAIL_CONTROL_EFFECT_NAME),
        PinKind::Starch => None,
        PinKind::Overlap => Some(OVERLAP_CONTROL_EFFECT_NAME),
    }
}

fn is_pin_control_effect(name: &str) -> bool {
    [
        POSITION_CONTROL_EFFECT_NAME,
        BONE_CONTROL_EFFECT_NAME,
        BEND_CONTROL_EFFECT_NAME,
        DETAIL_CONTROL_EFFECT_NAME,
        OVERLAP_CONTROL_EFFECT_NAME,
    ]
    .contains(&name)
}

fn control_has_xy(kind: PinKind) -> bool {
    matches!(kind, PinKind::Position | PinKind::Detail)
}

fn parse_control_pin_number(value: &str) -> Option<usize> {
    if let Ok(number) = value.trim().parse::<usize>() {
        return (number > 0).then_some(number);
    }
    let number = value.trim().parse::<f64>().ok()?;
    (number.is_finite() && number >= 1.0 && number.fract() == 0.0).then_some(number as usize)
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ControlUpdateAction {
    Add,
    AlreadyExists,
    Replace,
}

fn control_update_action(existing_names: &[String], desired_name: &str) -> ControlUpdateAction {
    match existing_names {
        [] => ControlUpdateAction::Add,
        [name] if name == desired_name => ControlUpdateAction::AlreadyExists,
        _ => ControlUpdateAction::Replace,
    }
}

fn add_pin_control(shared: &SharedState, target: EditTarget, pin_index: usize, kind: PinKind) {
    if shared.editor_snapshot().target.as_ref() != Some(&target) {
        return;
    }
    let Some(effect_name) = control_effect_name(kind) else {
        return;
    };
    shared.update_editor(|editor| {
        editor.busy = true;
        editor.status = format!("#{} {}制御を追加しています。", pin_index + 1, kind.label());
    });

    let target_for_edit = target.clone();
    let image_size = latest_image_size(shared);
    let result = EDIT_HANDLE.call_edit_section(move |edit| -> AnyResult<bool> {
        verify_target(edit.as_read_section(), &target_for_edit)?;
        let loaded = read_model(edit.as_read_section(), &target_for_edit, image_size)?;
        if loaded.model.kinds.get(pin_index) != Some(&kind) {
            aviutl2::anyhow::bail!("ピンの種類または番号が変更されています");
        }
        let source = loaded.model.src[pin_index];
        let mut existing_controls = Vec::new();
        let mut existing_names = Vec::new();
        for control in pin_control_effects_before_target(edit.as_read_section(), &target_for_edit)?
        {
            let value = edit
                .as_read_section()
                .get_effect_item_value(control, CONTROL_PIN_NUMBER_ITEM_NAME)
                .context("既存の制御フィルタ番号を取得できませんでした")?;
            if parse_control_pin_number(&value) == Some(pin_index + 1) {
                existing_controls.push(control);
                existing_names.push(
                    edit.as_read_section()
                        .get_effect_name(control)
                        .context("既存の制御フィルタ名を取得できませんでした")?,
                );
            }
        }
        let update_action = control_update_action(&existing_names, effect_name);
        if update_action == ControlUpdateAction::AlreadyExists {
            aviutl2::anyhow::bail!(
                "#{} の{}は既に追加されています",
                pin_index + 1,
                kind.label()
            );
        }
        let (_, target_position) =
            find_target_effect_position(edit.as_read_section(), &target_for_edit)?;
        let control = edit
            .create_effect(target_for_edit.object, effect_name)
            .with_context(|| format!("「{}」を追加できませんでした", kind.label()))?;
        let configure_result = (|| -> AnyResult<()> {
            let actual_position = edit
                .move_effect(target_for_edit.object, control, target_position)
                .context("制御フィルタをパペットツールの直前へ移動できませんでした")?;
            if actual_position != target_position {
                aviutl2::anyhow::bail!("制御フィルタをパペットツールの直前へ配置できませんでした");
            }
            edit.set_effect_item_value(
                control,
                CONTROL_PIN_NUMBER_ITEM_NAME,
                &(pin_index + 1).to_string(),
            )
            .context("制御フィルタへピン番号を書き込めませんでした")?;
            if control_has_xy(kind) {
                edit.set_effect_item_value(control, CONTROL_X_ITEM_NAME, &format_number(source.x))
                    .context("制御フィルタへピン元X座標を書き込めませんでした")?;
                edit.set_effect_item_value(control, CONTROL_Y_ITEM_NAME, &format_number(source.y))
                    .context("制御フィルタへピン元Y座標を書き込めませんでした")?;
            }
            for existing in &existing_controls {
                edit.delete_effect(target_for_edit.object, *existing)
                    .context("既存のピン制御フィルタを置き換えられませんでした")?;
            }
            Ok(())
        })();
        if let Err(error) = configure_result {
            let _ = edit.delete_effect(target_for_edit.object, control);
            return Err(error);
        }
        Ok(update_action == ControlUpdateAction::Replace)
    });

    match result {
        Ok(Ok(replaced)) => shared.update_editor(|editor| {
            if editor.target.as_ref() == Some(&target) {
                editor.busy = false;
                editor.status = if replaced {
                    format!("#{} {}制御を置き換えました。", pin_index + 1, kind.label())
                } else {
                    format!("#{} {}制御を追加しました。", pin_index + 1, kind.label())
                };
            }
        }),
        Ok(Err(error)) => pin_control_failed(shared, &target, format!("{error:#}")),
        Err(error) => pin_control_failed(shared, &target, error.to_string()),
    }
}

fn pin_control_failed(shared: &SharedState, target: &EditTarget, error: String) {
    shared.update_editor(|editor| {
        if editor.target.as_ref() == Some(target) {
            editor.busy = false;
            editor.status = format!("ピン制御を追加できませんでした: {error}");
        }
    });
}

fn add_object(shared: &SharedState) {
    shared.update_editor(|editor| {
        editor.busy = true;
        editor.status = "オブジェクトを追加しています。".to_owned();
    });

    let capture_after = shared.current_capture_sequence();
    let result = EDIT_HANDLE.call_edit_section(move |edit| -> AnyResult<EditTarget> {
        let object = edit.create_object(
            PUPPET_EFFECT_NAME,
            edit.info.layer,
            edit.info.frame,
            Some(NEW_OBJECT_LENGTH),
        )?;
        edit.set_focus_object(Some(object))?;
        if edit.count_object_effect(object, PUPPET_EFFECT_NAME)? == 0 {
            aviutl2::anyhow::bail!("追加したオブジェクトに対象効果がありません");
        }
        Ok(EditTarget {
            object,
            effect_name: PUPPET_EFFECT_NAME.to_owned(),
            effect_index: 0,
        })
    });

    match result {
        Ok(Ok(target)) => select_target(shared, target, capture_after),
        Ok(Err(error)) => shared.update_editor(|editor| {
            editor.busy = false;
            editor.status = format!("オブジェクトを追加できませんでした: {error:#}");
        }),
        Err(error) => shared.update_editor(|editor| {
            editor.busy = false;
            editor.status = format!("オブジェクトを追加できませんでした: {error}");
        }),
    }
}

fn select_focused(shared: &SharedState, effect_index: usize) {
    if !EDIT_HANDLE.is_ready() {
        return;
    }
    shared.update_editor(|editor| {
        editor.busy = true;
        editor.status = "選択オブジェクトを読み込んでいます。".to_owned();
    });
    let result = EDIT_HANDLE.call_read_section(move |read| -> AnyResult<EditTarget> {
        let object = read
            .get_focused_object()?
            .ok_or_else(|| aviutl2::anyhow::anyhow!("オブジェクトが選択されていません"))?;
        let effect_count = read.count_object_effect(object, PUPPET_EFFECT_NAME)?;
        if effect_count == 0 {
            aviutl2::anyhow::bail!("選択中のオブジェクトに対象効果がありません");
        }
        if effect_index >= effect_count {
            aviutl2::anyhow::bail!(
                "効果番号 {} は存在しません（対象効果は {} 個）",
                effect_index + 1,
                effect_count
            );
        }
        Ok(EditTarget {
            object,
            effect_name: PUPPET_EFFECT_NAME.to_owned(),
            effect_index,
        })
    });
    match result {
        Ok(Ok(target)) => {
            let capture_after = shared.current_capture_sequence();
            select_target(shared, target, capture_after);
        }
        Ok(Err(error)) => clear_target(shared, format!("{error:#}")),
        Err(error) => clear_target(shared, error.to_string()),
    }
}

fn clear_target(shared: &SharedState, status: String) {
    shared.update_editor(|editor| {
        editor.target = None;
        editor.model = None;
        editor.busy = false;
        editor.capture_after = shared.current_capture_sequence();
        editor.status = status;
    });
}

fn select_target(shared: &SharedState, target: EditTarget, capture_after: u64) {
    shared.update_editor(|editor| {
        editor.target = Some(target.clone());
        editor.model = None;
        editor.busy = true;
        editor.status = "対象のパラメーターを読み込んでいます。".to_owned();
        editor.capture_after = capture_after;
    });
    match load_model(shared, &target) {
        Ok(model) => shared.update_editor(|editor| {
            if editor.target.as_ref() == Some(&target) {
                editor.model = Some(model);
                editor.busy = false;
                editor.status = format!(
                    "{} 個のピンを編集中",
                    editor.model.as_ref().unwrap().src.len()
                );
            }
        }),
        Err(error) => shared.update_editor(|editor| {
            if editor.target.as_ref() == Some(&target) {
                editor.model = None;
                editor.busy = false;
                editor.status = format!("対象を開けません: {error:#}");
            }
        }),
    }
}

fn refresh(shared: &SharedState) {
    let Some(target) = shared.editor_snapshot().target else {
        select_focused(shared, 0);
        return;
    };
    match load_model(shared, &target) {
        Ok(model) => shared.update_editor(|editor| {
            if editor.target.as_ref() == Some(&target) && !editor.busy {
                editor.model = Some(model);
                editor.status = format!(
                    "{} 個のピンを編集中",
                    editor.model.as_ref().unwrap().src.len()
                );
            }
        }),
        Err(error) => shared.update_editor(|editor| {
            if editor.target.as_ref() == Some(&target) {
                editor.model = None;
                editor.busy = false;
                editor.status = format!("再読込に失敗しました: {error:#}");
            }
        }),
    }
}

fn apply(shared: &SharedState, target: EditTarget, operation: PinOperation) {
    if shared.editor_snapshot().target.as_ref() != Some(&target) {
        return;
    }
    shared.update_editor(|editor| {
        editor.busy = true;
        editor.status = "変更を保存しています。".to_owned();
    });

    let target_for_edit = target.clone();
    let image_size = latest_image_size(shared);
    let result = EDIT_HANDLE.call_edit_section(move |edit| -> AnyResult<PinModel> {
        verify_target(edit.as_read_section(), &target_for_edit)?;
        let loaded = read_model(edit.as_read_section(), &target_for_edit, image_size)?;
        let mut model = loaded.model;
        model
            .apply(operation)
            .map_err(aviutl2::anyhow::Error::msg)?;
        if let PinOperation::Delete(deleted_index) = operation {
            synchronize_pin_controls_after_delete(edit, &target_for_edit, deleted_index)?;
        }

        edit.set_object_effect_item(
            target_for_edit.object,
            &target_for_edit.effect_name,
            target_for_edit.effect_index,
            SOURCE_ITEM_NAME,
            &model.src_string(),
        )
        .context("「ピン元」の書き込みに失敗しました")?;
        edit.set_object_effect_item(
            target_for_edit.object,
            &target_for_edit.effect_name,
            target_for_edit.effect_index,
            DESTINATION_ITEM_NAME,
            &model.dst_string(),
        )
        .context("「ピン先」の書き込みに失敗しました")?;
        edit.set_object_effect_item(
            target_for_edit.object,
            &target_for_edit.effect_name,
            target_for_edit.effect_index,
            PIN_TYPES_ITEM_NAME,
            &model.pin_types_string(),
        )
        .context("「ピン種類」の書き込みに失敗しました")?;
        edit.set_object_effect_item(
            target_for_edit.object,
            &target_for_edit.effect_name,
            target_for_edit.effect_index,
            BONE_FOREST_ITEM_NAME,
            &model.bone_forest_string(),
        )
        .context("「ボーン親子関係」の書き込みに失敗しました")?;
        if !matches!(operation, PinOperation::MoveSourceAndDestination { .. })
            || loaded.synthesized_values
        {
            edit.set_object_effect_item(
                target_for_edit.object,
                &target_for_edit.effect_name,
                target_for_edit.effect_index,
                PINS_ITEM_NAME,
                &model.pin_count_string(),
            )
            .context("「ピン数」の書き込みに失敗しました")?;
        }
        Ok(model)
    });

    match result {
        Ok(Ok(model)) => shared.update_editor(|editor| {
            if editor.target.as_ref() == Some(&target) {
                editor.model = Some(model);
                editor.busy = false;
                editor.status = format!(
                    "{} 個のピンを編集中",
                    editor.model.as_ref().unwrap().src.len()
                );
            }
        }),
        Ok(Err(error)) => apply_failed(shared, &target, format!("{error:#}")),
        Err(error) => apply_failed(shared, &target, error.to_string()),
    }
}

fn apply_failed(shared: &SharedState, target: &EditTarget, error: String) {
    let reloaded = load_model(shared, target);
    shared.update_editor(|editor| {
        if editor.target.as_ref() == Some(target) {
            editor.busy = false;
            editor.model = reloaded.ok();
            editor.status = format!("変更を保存できませんでした: {error}");
        }
    });
}

fn load_model(shared: &SharedState, target: &EditTarget) -> AnyResult<PinModel> {
    if !EDIT_HANDLE.is_ready() {
        aviutl2::anyhow::bail!("編集APIが利用できません");
    }
    let target = target.clone();
    let image_size = latest_image_size(shared);
    EDIT_HANDLE.call_read_section(move |read| -> AnyResult<PinModel> {
        verify_target(read, &target)?;
        read_model(read, &target, image_size).map(|loaded| loaded.model)
    })?
}

fn latest_image_size(shared: &SharedState) -> Option<(usize, usize)> {
    let capture_after = shared.editor_snapshot().capture_after;
    shared
        .capture_snapshot()
        .latest
        .filter(|frame| frame.sequence > capture_after)
        .map(|frame| (frame.width, frame.height))
}

fn find_target_effect_position(
    read: &aviutl2::generic::ReadSection,
    target: &EditTarget,
) -> AnyResult<(aviutl2::generic::EffectHandle, usize)> {
    let effects = read
        .get_effects(target.object)
        .context("エフェクト一覧を取得できませんでした")?;
    let mut matching_index = 0;
    for (position, effect) in effects.into_iter().enumerate() {
        let name = read
            .get_effect_name(effect)
            .context("エフェクト名を取得できませんでした")?;
        if name == target.effect_name {
            if matching_index == target.effect_index {
                return Ok((effect, position));
            }
            matching_index += 1;
        }
    }
    aviutl2::anyhow::bail!("対象のパペットツールが見つかりません")
}

fn pin_control_effects_before_target(
    read: &aviutl2::generic::ReadSection,
    target: &EditTarget,
) -> AnyResult<Vec<aviutl2::generic::EffectHandle>> {
    let effects = read
        .get_effects(target.object)
        .context("エフェクト一覧を取得できませんでした")?;
    let (_, target_position) = find_target_effect_position(read, target)?;
    let mut controls = Vec::new();
    for &effect in effects[..target_position].iter().rev() {
        let name = read
            .get_effect_name(effect)
            .context("制御フィルタ名を取得できませんでした")?;
        if !is_pin_control_effect(&name) {
            break;
        }
        controls.push(effect);
    }
    Ok(controls)
}

fn synchronize_pin_controls_after_delete(
    edit: &aviutl2::generic::EditSection,
    target: &EditTarget,
    deleted_index: usize,
) -> AnyResult<()> {
    let deleted_number = deleted_index + 1;
    let controls = pin_control_effects_before_target(edit.as_read_section(), target)?;
    for control in controls {
        let value = edit
            .as_read_section()
            .get_effect_item_value(control, CONTROL_PIN_NUMBER_ITEM_NAME)
            .context("制御フィルタのピン番号を取得できませんでした")?;
        let Some(number) = parse_control_pin_number(&value) else {
            continue;
        };
        if number == deleted_number {
            edit.delete_effect(target.object, control)
                .context("削除したピンの制御フィルタを削除できませんでした")?;
        } else if number > deleted_number {
            edit.set_effect_item_value(
                control,
                CONTROL_PIN_NUMBER_ITEM_NAME,
                &(number - 1).to_string(),
            )
            .context("制御フィルタのピン番号を詰められませんでした")?;
        }
    }
    Ok(())
}

fn verify_target(read: &aviutl2::generic::ReadSection, target: &EditTarget) -> AnyResult<()> {
    if !read.object_exists(target.object) {
        aviutl2::anyhow::bail!("オブジェクトが存在しません");
    }
    let effect_count = read
        .count_object_effect(target.object, &target.effect_name)
        .context("対象効果の確認に失敗しました")?;
    if target.effect_index >= effect_count {
        aviutl2::anyhow::bail!("対象効果が削除されたか、効果番号が変わりました");
    }
    Ok(())
}

struct LoadedModel {
    model: PinModel,
    synthesized_values: bool,
}

fn read_model(
    read: &aviutl2::generic::ReadSection,
    target: &EditTarget,
    image_size: Option<(usize, usize)>,
) -> AnyResult<LoadedModel> {
    let pins = read
        .get_object_effect_item(
            target.object,
            &target.effect_name,
            target.effect_index,
            PINS_ITEM_NAME,
        )
        .context("「ピン数」の読み込みに失敗しました")?;
    let src = read
        .get_object_effect_item(
            target.object,
            &target.effect_name,
            target.effect_index,
            SOURCE_ITEM_NAME,
        )
        .context("「ピン元」の読み込みに失敗しました")?;
    let dst = read
        .get_object_effect_item(
            target.object,
            &target.effect_name,
            target.effect_index,
            DESTINATION_ITEM_NAME,
        )
        .context("「ピン先」の読み込みに失敗しました")?;
    let pin_types = read
        .get_object_effect_item(
            target.object,
            &target.effect_name,
            target.effect_index,
            PIN_TYPES_ITEM_NAME,
        )
        .context("「ピン種類」の読み込みに失敗しました")?;
    let bone_forest = read
        .get_object_effect_item(
            target.object,
            &target.effect_name,
            target.effect_index,
            BONE_FOREST_ITEM_NAME,
        )
        .context("「ボーン親子関係」の読み込みに失敗しました")?;
    let (model, synthesized_values) =
        PinModel::parse_or_default(&pins, &src, &dst, &pin_types, &bone_forest, image_size)
            .map_err(aviutl2::anyhow::Error::msg)?;
    Ok(LoadedModel {
        model,
        synthesized_values,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_supported_pin_kind_maps_to_its_control_effect() {
        assert_eq!(
            control_effect_name(PinKind::Position),
            Some(POSITION_CONTROL_EFFECT_NAME)
        );
        assert_eq!(
            control_effect_name(PinKind::Bone),
            Some(BONE_CONTROL_EFFECT_NAME)
        );
        assert_eq!(
            control_effect_name(PinKind::Bend),
            Some(BEND_CONTROL_EFFECT_NAME)
        );
        assert_eq!(
            control_effect_name(PinKind::Detail),
            Some(DETAIL_CONTROL_EFFECT_NAME)
        );
        assert_eq!(control_effect_name(PinKind::Starch), None);
        assert_eq!(
            control_effect_name(PinKind::Overlap),
            Some(OVERLAP_CONTROL_EFFECT_NAME)
        );
    }

    #[test]
    fn only_known_control_effects_are_grouped_with_a_puppet_effect() {
        for kind in PinKind::ALL {
            if let Some(name) = control_effect_name(kind) {
                assert!(is_pin_control_effect(name));
            }
        }
        assert!(!is_pin_control_effect(PUPPET_EFFECT_NAME));
        assert!(!is_pin_control_effect("ぼかし"));
    }

    #[test]
    fn xy_is_initialized_only_for_controls_that_expose_xy_tracks() {
        assert!(control_has_xy(PinKind::Position));
        assert!(!control_has_xy(PinKind::Bone));
        assert!(!control_has_xy(PinKind::Bend));
        assert!(control_has_xy(PinKind::Detail));
        assert!(!control_has_xy(PinKind::Starch));
        assert!(!control_has_xy(PinKind::Overlap));
    }

    #[test]
    fn control_pin_numbers_accept_integer_track_representations() {
        assert_eq!(parse_control_pin_number("3"), Some(3));
        assert_eq!(parse_control_pin_number("3.0"), Some(3));
        assert_eq!(parse_control_pin_number("0"), None);
        assert_eq!(parse_control_pin_number("3.5"), None);
        assert_eq!(parse_control_pin_number("invalid"), None);
    }

    #[test]
    fn control_updates_compare_both_pin_id_and_effect_type() {
        let desired = POSITION_CONTROL_EFFECT_NAME;
        assert_eq!(
            control_update_action(&[], desired),
            ControlUpdateAction::Add
        );
        assert_eq!(
            control_update_action(&[desired.to_owned()], desired),
            ControlUpdateAction::AlreadyExists
        );
        assert_eq!(
            control_update_action(&[BONE_CONTROL_EFFECT_NAME.to_owned()], desired),
            ControlUpdateAction::Replace
        );
        assert_eq!(
            control_update_action(&[desired.to_owned(), desired.to_owned()], desired),
            ControlUpdateAction::Replace
        );
    }
}
