-- 測地距離・MLS・細分化は、一つの変形パイプラインとしてRust側で処理する。
local mesh_vertices, mesh_indices = {}, {}
for i = 1, mesh_n_verts do
	local n = #mesh_vertices
	mesh_vertices[n + 1], mesh_vertices[n + 2] = mesh_x[i], mesh_y[i]
end
for _, tri in ipairs(mesh_tris) do
	for j = 1, 3 do mesh_indices[#mesh_indices + 1] = tri[j] - 1 end
end
local pin_data = {}
local pin_dynamics_debug = {}
for i = 1, pins do
	pin_data[i] = {
		sx = pin_sx[i], sy = pin_sy[i], dx = pin_dx[i], dy = pin_dy[i],
		kind = pin_types[i], layer = pin_layer[i],
		rotation = pin_rotation[i], scale = pin_scale[i],
		range = pin_range[i], show_range = pin_show_range[i],
	}
end
---$include "PuppetPinDynamics.inc.lua"

local function run_deformation(sources, destinations, layers, ranges, force_mls)
	-- モジュール取得失敗も下のpcallで簡易変形へフォールバックさせる。
	local deform = not force_mls and deformationMethod == 2
		and mesh_module.deform_arap or mesh_module.deform_mls
	local layer_payload = {}
	for i = 1, #layers do layer_payload[i] = layers[i] end
	for i = 1, #ranges do layer_payload[#layers + i] = ranges[i] end
	return deform(mesh_vertices, mesh_indices, sources, destinations,
		layer_payload, stiff, div, w, h)
end


local ok, deformed, render_vertices = pcall(function()
	-- まず移動先を持つピンだけで自然な姿勢を作る。ベンド/スターチを
	-- 元位置のまま混ぜるとアンカーになり、追従できなくなるため分離する。
	local base_sources, base_destinations, base_layers, base_ranges = {}, {}, {}, {}
	for _, pin in ipairs(pin_data) do
		if pin.kind == PIN_TYPE.POSITION or pin.kind == PIN_TYPE.BONE
			or pin.kind == PIN_TYPE.DETAIL then
			local n = #base_sources
			base_sources[n + 1], base_sources[n + 2] = pin.sx, pin.sy
			base_destinations[n + 1], base_destinations[n + 2] = pin.dx, pin.dy
			base_layers[#base_layers + 1] = pin.layer
			base_ranges[#base_ranges + 1] = 0
		end
	end
	if #base_sources == 0 then
		-- 0ピン時は先頭メッシュ頂点を動かさない内部拘束として使う。
		-- UI上のピンではなく、変形モジュールの非空入力要件を満たすだけ。
		local anchor_x, anchor_y = mesh_vertices[1], mesh_vertices[2]
		base_sources = { anchor_x, anchor_y }
		base_destinations = { anchor_x, anchor_y }
		base_layers = { 0 }
		base_ranges = { 0 }
	end
	local needs_dynamic_pose = false
	for _, pin in ipairs(pin_data) do
		if pin.kind == PIN_TYPE.BEND or pin.kind == PIN_TYPE.DETAIL
			or pin.kind == PIN_TYPE.STARCH or pin.kind == PIN_TYPE.OVERLAP then
			needs_dynamic_pose = true
			break
		end
	end
	if not needs_dynamic_pose then
		return run_deformation(base_sources, base_destinations, base_layers, base_ranges)
	end
	local preliminary = run_deformation(
		base_sources, base_destinations, base_layers, base_ranges)
	-- これは最小値だけ。実際の影響範囲は各ピンに最も近いメッシュ領域の
	-- 外周まで自動的に拡張される。
	local handle_radius = math.max(6, math.max(w, h) / math.max(density, 1))
	local pin_sources, pin_destinations, final_layers, final_ranges, debug_handles =
		puppet_pin_dynamics.build(
		mesh_vertices, preliminary, pin_data, handle_radius, PIN_TYPE)
	pin_dynamics_debug = debug_handles
	if #pin_sources == 0 then
		pin_sources, pin_destinations, final_layers, final_ranges =
			base_sources, base_destinations, base_layers, base_ranges
	end
	local has_geometry_pin = false
	for _, encoded_range in ipairs(final_ranges) do
		if encoded_range >= 0 then has_geometry_pin = true; break end
	end
	if not has_geometry_pin then
		local n = #pin_sources
		pin_sources[n + 1], pin_sources[n + 2] = base_sources[1], base_sources[2]
		pin_destinations[n + 1], pin_destinations[n + 2] =
			base_destinations[1], base_destinations[2]
		final_layers[#final_layers + 1] = 0
		final_ranges[#final_ranges + 1] = 0
	end
	local final_deformed, final_render = run_deformation(
		pin_sources, pin_destinations, final_layers, final_ranges)
	return final_deformed, final_render
end)
if not ok then
	pin_dynamics_debug = {}
	if not puppet_deformation_error_reported then
		print("@warn", "puppet_geometry: 変形に失敗: " .. tostring(deformed) .. " (簡易変形を使用します)")
		puppet_deformation_error_reported = true
	end
	deformed, render_vertices = {}, {}
	for i = 1, mesh_n_verts do
		local x, y = mesh_x[i], mesh_y[i]
		if pins == 0 then
			deformed[i*2-1], deformed[i*2] = x, y
		else
			local dx, dy, total = 0, 0, 0
			for p = 1, pins do
				local distance = math.max(math.sqrt((x-pin_sx[p])^2 + (y-pin_sy[p])^2), 1e-5)
				local weight = distance ^ (-2 * stiff)
				dx, dy, total = dx + weight*(pin_dx[p]-pin_sx[p]), dy + weight*(pin_dy[p]-pin_sy[p]), total + weight
			end
			deformed[i*2-1], deformed[i*2] = x + dx/total, y + dy/total
		end
	end
	for _, tri in ipairs(mesh_tris) do
		for j = 1, 3 do
			local i = tri[j]
			local n = #render_vertices
			render_vertices[n+1], render_vertices[n+2] = deformed[i*2-1], deformed[i*2]
			render_vertices[n+3], render_vertices[n+4] = (mesh_x[i]+hw)/w, (mesh_y[i]+hh)/h
		end
	end
else
	puppet_deformation_error_reported = false
end
local def_x, def_y = {}, {}
for i = 1, #deformed, 2 do
	def_x[(i + 1) / 2], def_y[(i + 1) / 2] = deformed[i], deformed[i + 1]
end
local show_gui = show and obj.getoption("gui")
