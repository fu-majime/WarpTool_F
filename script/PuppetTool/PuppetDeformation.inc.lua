-- MLS/ARAPの評価と測地距離と描画処理はRust側で行う
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
local handle_radius = math.max(1, math.max(w, h) * 0.015)
local material_pins = {}
---$include "PuppetPinDynamics.inc.lua"
for i = 1, pins do
	pin_data[i] = {
		sx = pin_sx[i], sy = pin_sy[i], dx = pin_dx[i], dy = pin_dy[i],
		kind = pin_types[i], layer = pin_layer[i],
		rotation = pin_rotation[i], scale = pin_scale[i],
		range = pin_range[i], show_range = pin_show_range[i],
	}
end
local bend_frames = puppet_pin_dynamics.bend_frames(pin_data, PIN_TYPE)
for i = 1, pins do
	if pin_types[i] ~= PIN_TYPE.OVERLAP and
		(pin_types[i] ~= PIN_TYPE.BEND or bend_frames or math.abs(pin_rotation[i]) > 1e-8 or math.abs(pin_scale[i]-1) > 1e-8) then
		material_pins[#material_pins + 1] = pin_sx[i]
		material_pins[#material_pins + 1] = pin_sy[i]
	end
end
if #mesh_indices > 0 then
	mesh_vertices, mesh_indices = mesh_module.prepare_mesh(
		mesh_vertices, mesh_indices, material_pins, handle_radius, math.max(w,h)/density)
	mesh_x, mesh_y, mesh_tris = {}, {}, {}
	mesh_n_verts = #mesh_vertices / 2
	for i = 1, mesh_n_verts do
		mesh_x[i], mesh_y[i] = mesh_vertices[i*2-1], mesh_vertices[i*2]
	end
	for i = 1, #mesh_indices, 3 do
		mesh_tris[#mesh_tris + 1] = { mesh_indices[i]+1, mesh_indices[i+1]+1, mesh_indices[i+2]+1 }
	end
end

local function run_deformation(sources, destinations, layers, ranges, vertices_only, poses)
	local deform
	if deformationMethod == 2 then deform = mesh_module.deform_arap
	else deform = mesh_module.deform_mls end
	local layer_payload = {}
	for i = 1, #layers do layer_payload[i] = layers[i] end
	for i = 1, #ranges do layer_payload[#layers + i] = ranges[i] end
	local render_divisions = vertices_only and 0 or (deformationMethod == 2 and 1 or div)
	return deform(mesh_vertices, mesh_indices, sources, destinations,
		layer_payload, stiff, render_divisions, w, h, poses or {})
end


local ok, deformed, render_vertices, wire_vertices = pcall(function()
	if #mesh_indices == 0 then return {}, {}, {} end
	-- まず移動先を持つピンだけで自然な姿勢を作る。ベンド/スターチは混ぜないらしい
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
	local needs_dynamic_pose = #base_sources > 2 or pins > 0
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
	-- ARAP derives position/bone rotations in the native solve, and
	-- detail controls prescribe their angle. Only following controls need
	-- a preliminary material position; avoid a second full solve otherwise.
	local needs_preliminary = deformationMethod ~= 2
	for _, pin in ipairs(pin_data) do
		if pin.kind == PIN_TYPE.BEND or pin.kind == PIN_TYPE.STARCH
			or pin.kind == PIN_TYPE.OVERLAP then needs_preliminary = true; break end
	end
	local preliminary = mesh_vertices
	if needs_preliminary then
		preliminary = run_deformation(
			base_sources, base_destinations, base_layers, base_ranges, true)
	end
	local pin_sources, pin_destinations, final_layers, final_ranges, debug_handles, poses =
		puppet_pin_dynamics.build(
		mesh_vertices, preliminary, pin_data, handle_radius, PIN_TYPE, mesh_indices)
	pin_dynamics_debug = debug_handles
	local final_deformed, final_render, final_wire = run_deformation(
		pin_sources, pin_destinations, final_layers, final_ranges, false, poses)
	return final_deformed, final_render, final_wire
end)
if not ok then
	error("puppet_geometry: 選択した変形方式に失敗: " .. tostring(deformed))
end
local def_x, def_y = {}, {}
for i = 1, #deformed, 2 do
	def_x[(i + 1) / 2], def_y[(i + 1) / 2] = deformed[i], deformed[i + 1]
end
local is_gui = obj.getoption("gui")
local show_gui = show and is_gui
local show_overlap_range = false
if is_gui then
	for _, pin in ipairs(pin_data) do
		if pin.kind == PIN_TYPE.OVERLAP and pin.show_range and pin.range > 0 then
			show_overlap_range = true
			break
		end
	end
end
if show_gui then
	for _,handle in ipairs(pin_dynamics_debug) do
		local u,v=handle.sx/w+0.5,handle.sy/h+0.5
		for i=1,#render_vertices,12 do
			local ax,ay=render_vertices[i+2],render_vertices[i+3]
			local bx,by=render_vertices[i+6]-ax,render_vertices[i+7]-ay
			local cx,cy=render_vertices[i+10]-ax,render_vertices[i+11]-ay
			local det=bx*cy-by*cx
			if math.abs(det)>1e-15 then
				local b=((u-ax)*cy-(v-ay)*cx)/det
				local c=(bx*(v-ay)-by*(u-ax))/det
				local a=1-b-c
				if a>=-1e-7 and b>=-1e-7 and c>=-1e-7 then
					handle.dx=a*render_vertices[i]+b*render_vertices[i+4]+c*render_vertices[i+8]
					handle.dy=a*render_vertices[i+1]+b*render_vertices[i+5]+c*render_vertices[i+9]
					break
				end
			end
		end
	end
end
