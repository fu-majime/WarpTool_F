-- =============================================================
-- 1. ピンのデフォルト位置（画像上に均等配置）
--    ユーザーがドラッグで関節等へ再配置する
-- =============================================================
local function make_default_pin_pos(n)
	local t = {}
	if n <= 0 then
		return t
	elseif n == 1 then
		t[1] = 0; t[2] = 0
	elseif n == 2 then
		t[1] = -hw * 0.3; t[2] = 0
		t[3] =  hw * 0.3; t[4] = 0
	else
		-- 中央1点 + 残りを楕円上に分散配置
		t[1] = 0; t[2] = 0
		local rx = hw * 0.4
		local ry = hh * 0.4
		for i = 1, n - 1 do
			local angle = (i - 1) / (n - 1) * math.pi * 2 - math.pi / 2
			t[i * 2 + 1] = math.cos(angle) * rx
			t[i * 2 + 2] = math.sin(angle) * ry
		end
	end
	return t
end

local default_pins = make_default_pin_pos(pins)

-- 移動先を持つピンだけを、元のピンID順にピン先へ格納する。
-- ベンド・スターチ・重なりピンは移動先を持たない。
local destination_pin_ids = {}
local default_destinations = {}
for pin_id = 1, pins do
	local pin_type = pin_types[pin_id]
	if pin_type ~= PIN_TYPE.BEND
		and pin_type ~= PIN_TYPE.STARCH
		and pin_type ~= PIN_TYPE.OVERLAP then
		destination_pin_ids[#destination_pin_ids + 1] = pin_id
		default_destinations[#default_destinations + 1] = default_pins[pin_id * 2 - 1]
		default_destinations[#default_destinations + 1] = default_pins[pin_id * 2]
	end
end

-- ピン数変更時はリセット
if #src ~= #default_pins then src = default_pins end
default_destinations = {}
for _, pin_id in ipairs(destination_pin_ids) do
	default_destinations[#default_destinations + 1] = src[pin_id * 2 - 1]
	default_destinations[#default_destinations + 1] = src[pin_id * 2]
end
if #dst == #default_pins and #dst ~= #default_destinations then
	-- 旧形式（全ピン分の移動先）をID順を保ったまま圧縮する。
	local legacy_dst = dst
	dst = {}
	for _, pin_id in ipairs(destination_pin_ids) do
		dst[#dst + 1] = legacy_dst[pin_id * 2 - 1]
		dst[#dst + 1] = legacy_dst[pin_id * 2]
	end
elseif #dst ~= #default_destinations then
	-- dstのデフォルトは、移動先を持つ各srcと同じ位置（変形なし）
	dst = {}
	for _, pin_id in ipairs(destination_pin_ids) do
		dst[#dst + 1] = src[pin_id * 2 - 1]
		dst[#dst + 1] = src[pin_id * 2]
	end
end

-- アンカー設定: チェック中は変形元、チェックを外すと変形先を操作する TODO ピン制御時は使わないときもあるはず
if show and pins > 0 then
    if edit_source then
        obj.setanchor("src", pins, default_pins, "rgba", 0x20e060e0)
    else
        obj.setanchor("dst", #destination_pin_ids, default_destinations, "rgba", 0xf05050e0)
    end
end

-- 変形率を反映
local pin_sx, pin_sy = {}, {}  -- 変形元
local pin_tx, pin_ty = {}, {}  -- 変形先（ratio反映前）
local pin_dx, pin_dy = {}, {}  -- 変形先（ratio反映済み）
local pin_layer = {}           -- ユーザーが後から整数を足し算してレイヤー調整できるよう用意
local pin_rotation = {}        -- ベンド/詳細ピンの追加回転（degree）
local pin_scale = {}           -- ベンド/詳細ピンの追加拡大率（1.0 = 等倍）
local pin_range = {}           -- 重なりピンの影響半径
local pin_show_range = {}      -- 重なり範囲の可視化

-- ピン制御から値を取得
local effect_overrides = {}
local bone_overrides = {}
local bend_overrides = {}
local overlap_overrides = {}

local function collect_pin_effects(effect_label, read_position, read_bend)
	local fx_idx = 0
	while true do
		local effect_name = effect_label .. "@${PACKAGE_NAME}:" .. fx_idx
		local ok, pin_num = pcall(function() return obj.getvalue(effect_name, "番号") end)
		if not ok or pin_num == nil then break end
		local pin_id = math.floor(tonumber(pin_num) or 0)
		if pin_id >= 1 and pin_id <= pins then
			if read_position then
				local ok_x, pin_x = pcall(function() return obj.getvalue(effect_name, "X") end)
				local ok_y, pin_y = pcall(function() return obj.getvalue(effect_name, "Y") end)
				if ok_x and ok_y and tonumber(pin_x) and tonumber(pin_y) then
					effect_overrides[pin_id] = { tonumber(pin_x), tonumber(pin_y) }
				end
			end
			if read_bend then
				local ok_rot, rot = pcall(function() return obj.getvalue(effect_name, "回転") end)
				local ok_scale, scale = pcall(function() return obj.getvalue(effect_name, "拡大率") end)
				if ok_rot and ok_scale and tonumber(rot) and tonumber(scale) then
					bend_overrides[pin_id] = { tonumber(rot), tonumber(scale) }
				end
			end
		end
		fx_idx = fx_idx + 1
	end
end

collect_pin_effects("${EFFECT_NAME_position_pin}", true, false)
collect_pin_effects("${EFFECT_NAME_bend_pin}", false, true)
collect_pin_effects("${EFFECT_NAME_advanced_pin}", true, true)

local bone_fx_idx = 0
while true do
	local effect_name = "${EFFECT_NAME_bone_pin}@${PACKAGE_NAME}:" .. bone_fx_idx
	local ok, pin_num = pcall(function() return obj.getvalue(effect_name, "番号") end)
	if not ok or pin_num == nil then break end
	local pin_id = math.floor(tonumber(pin_num) or 0)
	if pin_id >= 1 and pin_id <= pins and pin_types[pin_id] == PIN_TYPE.BONE then
		local ok_rot, rotation = pcall(function() return obj.getvalue(effect_name, "回転") end)
		local ok_stretch, stretch = pcall(function() return obj.getvalue(effect_name, "伸縮") end)
		local ok_apply, apply = pcall(function() return obj.getvalue(effect_name, "ピン先に反映") end)
		if ok_rot and ok_stretch and tonumber(rotation) and tonumber(stretch) then
			bone_overrides[pin_id] = {
				rotation = tonumber(rotation), stretch = tonumber(stretch) / 100,
				apply = ok_apply and (apply == true or tonumber(apply) == 1),
			}
		end
	end
	bone_fx_idx = bone_fx_idx + 1
end

local fx_idx = 0
while true do
	local effect_name = "${EFFECT_NAME_overlap_pin}@${PACKAGE_NAME}:" .. fx_idx
	local ok, pin_num = pcall(function() return obj.getvalue(effect_name, "番号") end)
	if not ok or pin_num == nil then break end
	local pin_id = math.floor(tonumber(pin_num) or 0)
	if pin_id >= 1 and pin_id <= pins then
		local ok_front, front = pcall(function() return obj.getvalue(effect_name, "前方[%]") end)
		local ok_range, range_value = pcall(function() return obj.getvalue(effect_name, "範囲") end)
		local ok_show, show_range = pcall(function() return obj.getvalue(effect_name, "範囲表示") end)
		if ok_front and ok_range and tonumber(front) and tonumber(range_value) then
			overlap_overrides[pin_id] = {
				front = tonumber(front), range = math.max(0, tonumber(range_value)),
				show = ok_show and (show_range == true or tonumber(show_range) == 1),
			}
		end
	end
	fx_idx = fx_idx + 1
end

-- アンカーで指定した姿勢を解いた後、存在するボーン制御を優先して上書きする。
local bone_anchor_x, bone_anchor_y = {}, {}
local bone_target_x, bone_target_y, bone_write_control = {}, {}, {}
local function parse_number_list(value)
	local out = {}
	if type(value) == "string" then
		for token in value:gmatch("[^,]+") do out[#out + 1] = tonumber(token) end
	end
	return out
end
local function encode_number_list(values)
	local out = {}
	for i = 1, #values do out[i] = string.format("%.17g", values[i]) end
	return table.concat(out, ",")
end
local function encode_bone_controls(values)
	local out = {}
	for pin_id = 1, pins do
		local control = values[pin_id]
		out[#out + 1] = control and control.rotation or 0
		out[#out + 1] = control and control.stretch or 1
	end
	return encode_number_list(out)
end
local function decode_bone_controls(value)
	local numbers, out = parse_number_list(value), {}
	for pin_id = 1, pins do
		out[pin_id] = {
			rotation = numbers[pin_id * 2 - 1] or 0,
			stretch = numbers[pin_id * 2] or 1,
		}
	end
	return out
end
if not edit_source and type(bone_forest) == "table" then
	local forest = bone_forest
	while #forest == 1 and type(forest[1]) == "table"
		and type(forest[1][1]) == "table" do forest = forest[1] end
	local full_dst = {}
	for pin_id = 1, pins do
		full_dst[pin_id * 2 - 1], full_dst[pin_id * 2] = src[pin_id * 2 - 1], src[pin_id * 2]
	end
	for slot, pin_id in ipairs(destination_pin_ids) do
		full_dst[pin_id * 2 - 1], full_dst[pin_id * 2] = dst[slot * 2 - 1], dst[slot * 2]
	end
	-- 制御付きピンはアンカーを掴んだ瞬間だけ生じる座標も入力にしない。
	-- 最後に確定したピン先を使い、子ボーンの局所角度への混入を防ぐ。
	local last_controlled_dst = parse_number_list(global.puppet_bone_control_last)
	if #last_controlled_dst == #dst then
		for slot, pin_id in ipairs(destination_pin_ids) do
			if bone_overrides[pin_id] then
				full_dst[pin_id * 2 - 1] = last_controlled_dst[slot * 2 - 1]
				full_dst[pin_id * 2] = last_controlled_dst[slot * 2]
			end
		end
	end
	local applied_controls = decode_bone_controls(global.puppet_bone_control_applied)
	local pending = parse_number_list(global.puppet_bone_control_pending)
	if #pending > 0 then
		local arrived = #pending == #dst
		for i = 1, #dst do
			if math.abs((tonumber(dst[i]) or 0) - pending[i]) >= 0.1 then
				arrived = false
				break
			end
		end
		if arrived then
			global.puppet_bone_control_applied = global.puppet_bone_control_pending_controls
			applied_controls = decode_bone_controls(global.puppet_bone_control_applied)
			global.puppet_bone_control_pending = ""
			global.puppet_bone_control_pending_controls = ""
		else
			-- Bridge反映前の古いdstへ同じ制御差分を再適用しない。
			-- 待機中は要求済みの姿勢全体を次の計算入力として固定する。
			applied_controls = decode_bone_controls(global.puppet_bone_control_pending_controls)
			for slot, pin_id in ipairs(destination_pin_ids) do
				full_dst[pin_id * 2 - 1] = pending[slot * 2 - 1]
				full_dst[pin_id * 2] = pending[slot * 2]
			end
		end
	end
	local has_bone_overrides = next(bone_overrides) ~= nil
	local solved, state = pin_hierarchy.solve(
		src, full_dst, forest, pins,
		has_bone_overrides and nil or global.puppet_bone_hierarchy_state, false)
	if has_bone_overrides then
		-- 制御後の保存座標はアンカー解決座標と一致しないため、
		-- ソルバーの非同期待機フラグを持ち越さない。
		global.puppet_bone_hierarchy_state = ""
	else
		global.puppet_bone_hierarchy_state = state
	end

	local function previous_control(pin_id)
		return applied_controls[pin_id]
	end
	local function pose(tree, parent_id, absolute_rotation, rotation_delta,
		root_stretch_delta, inherited_write)
		if type(tree) ~= "table" then return end
		local pin_id = math.floor(tonumber(tree[1]) or 0)
		if pin_id < 1 or pin_id > pins then return end
		bone_anchor_x[pin_id], bone_anchor_y[pin_id] = solved[pin_id * 2 - 1], solved[pin_id * 2]
		local control = bone_overrides[pin_id]
		inherited_write = inherited_write or (control and control.apply) or false
		bone_write_control[pin_id] = inherited_write
		if not parent_id then
			bone_target_x[pin_id], bone_target_y[pin_id] = solved[pin_id * 2 - 1], solved[pin_id * 2]
			if control then
				local previous = previous_control(pin_id)
				absolute_rotation = math.rad(control.rotation)
				rotation_delta = math.rad(control.rotation - previous.rotation)
				root_stretch_delta = control.stretch / math.max(previous.stretch, 1e-6)
			else
				absolute_rotation, rotation_delta, root_stretch_delta = 0, 0, nil
			end
		else
			local manual_dx = solved[pin_id * 2 - 1] - solved[parent_id * 2 - 1]
			local manual_dy = solved[pin_id * 2] - solved[parent_id * 2]
			local angle, length
			if control then
				local previous = previous_control(pin_id)
				local source_dx = src[pin_id * 2 - 1] - src[parent_id * 2 - 1]
				local source_dy = src[pin_id * 2] - src[parent_id * 2]
				absolute_rotation = absolute_rotation + math.rad(control.rotation)
				rotation_delta = rotation_delta + math.rad(control.rotation - previous.rotation)
				angle = math.atan2(source_dy, source_dx) + absolute_rotation
				length = math.sqrt(source_dx * source_dx + source_dy * source_dy) * control.stretch
			else
				angle = math.atan2(manual_dy, manual_dx) + rotation_delta
				length = math.sqrt(manual_dx * manual_dx + manual_dy * manual_dy)
			end
			if root_stretch_delta then length = length * root_stretch_delta end
			bone_target_x[pin_id] = bone_target_x[parent_id] + math.cos(angle) * length
			bone_target_y[pin_id] = bone_target_y[parent_id] + math.sin(angle) * length
			root_stretch_delta = nil
		end
		for i = 2, #tree do
			pose(tree[i], pin_id, absolute_rotation, rotation_delta,
				root_stretch_delta, inherited_write)
		end
	end
	for i = 1, #forest do pose(forest[i], nil, 0, 0, nil, false) end
end

-- アンカー操作とボーン制御を反映した最終位置を「ピン先」へ保存する。
if obj.getoption("gui") then
	local solved_dst, changed = {}, false
	for i = 1, #dst do solved_dst[i] = dst[i] end
	for slot, pin_id in ipairs(destination_pin_ids) do
		local x, y
		if bone_write_control[pin_id] then
			x, y = bone_target_x[pin_id], bone_target_y[pin_id]
		else
			x, y = bone_anchor_x[pin_id], bone_anchor_y[pin_id]
		end
		if x and y then
			local x_index = slot * 2 - 1
			if math.abs(solved_dst[x_index] - x) >= 0.1
				or math.abs(solved_dst[x_index + 1] - y) >= 0.1 then
				solved_dst[x_index], solved_dst[x_index + 1] = x, y
				changed = true
			end
		end
	end
	if changed then
		global.puppet_bone_control_last = encode_number_list(solved_dst)
		global.puppet_bone_control_pending = encode_number_list(solved_dst)
		global.puppet_bone_control_pending_controls = encode_bone_controls(bone_overrides)
		local ok, bridge = pcall(function() return obj.module("ScriptEditBridge") end)
		if ok and type(bridge) == "table" and type(bridge.request) == "function" then
			pcall(bridge.request, obj.getoption("script_name"), 0,
				"ピン先", pin_hierarchy.encode(solved_dst))
		end
	elseif type(global.puppet_bone_control_pending) ~= "string" then
		global.puppet_bone_control_last = encode_number_list(solved_dst)
		global.puppet_bone_control_applied = encode_bone_controls(bone_overrides)
	end
end

local destination_index = 0
for i = 1, pins do
	pin_layer[i] = 0
	pin_rotation[i] = 0
	pin_scale[i] = 1
	pin_range[i] = 0
	pin_show_range[i] = false
	pin_sx[i] = src[i * 2 - 1]
	pin_sy[i] = src[i * 2]

	local base_dx = pin_sx[i]
	local base_dy = pin_sy[i]
	local pin_type = pin_types[i]
	if pin_type ~= PIN_TYPE.BEND
		and pin_type ~= PIN_TYPE.STARCH
		and pin_type ~= PIN_TYPE.OVERLAP then
		destination_index = destination_index + 1
		base_dx = dst[destination_index * 2 - 1]
		base_dy = dst[destination_index * 2]
	end

	if effect_overrides[i] then
		base_dx = effect_overrides[i][1]
		base_dy = effect_overrides[i][2]
	end
	if bone_target_x[i] then
		base_dx, base_dy = bone_target_x[i], bone_target_y[i]
	end
	if bend_overrides[i] then
		pin_rotation[i] = bend_overrides[i][1] * ratio / 100
		pin_scale[i] = 1 + (bend_overrides[i][2] / 100 - 1) * ratio / 100
	end
	if overlap_overrides[i] then
		pin_layer[i] = overlap_overrides[i].front / 100 * ratio / 100
		pin_range[i] = overlap_overrides[i].range
		pin_show_range[i] = overlap_overrides[i].show
	end
	pin_tx[i] = base_dx
	pin_ty[i] = base_dy

	if ratio < 100 then
		pin_dx[i] = src[i * 2 - 1] + (pin_tx[i] - src[i * 2 - 1]) * ratio / 100
		pin_dy[i] = src[i * 2]     + (pin_ty[i] - src[i * 2])     * ratio / 100
	else
		pin_dx[i] = pin_tx[i]
		pin_dy[i] = pin_ty[i]
	end
end

-- =============================================================
-- ユーティリティ
-- =============================================================
local function check_alpha_at(px_x, px_y)
	local px = math.floor(px_x + hw)
	local py = math.floor(px_y + hh)
	px = math.max(0, math.min(w - 1, px))
	py = math.max(0, math.min(h - 1, py))
	local _, a = obj.getpixel(px, py, "col")
	return a * 255 >= threshold
end
