-- =============================================================
-- 2. メッシュ生成
-- =============================================================
-- mesh_x[], mesh_y[] = 全頂点座標
-- mesh_tris[] = {{ i1, i2, i3 }, ...}  (1-based index into mesh_x/y)
local mesh_x, mesh_y = {}, {}
local mesh_tris = {}
local mesh_n_verts = 0

-- AviUtl2スクリプトモジュールからLua配列としてメッシュを受け取る。
-- FFIポインタやRust側メモリの所有権をスクリプトへ持ち出さない。
local function import_module_mesh(vertices, indices)
	if type(vertices) ~= "table" or type(indices) ~= "table" or
		#vertices < 6 or #vertices % 2 ~= 0 or
		#indices < 3 or #indices % 3 ~= 0 then
		return false
	end

	local vertex_count = #vertices / 2
	for i = 1, #vertices, 2 do
		local x, y = vertices[i], vertices[i + 1]
		if type(x) ~= "number" or type(y) ~= "number" then return false end
		mesh_x[(i + 1) / 2] = x
		mesh_y[(i + 1) / 2] = y
	end

	for i = 1, #indices, 3 do
		local i0, i1, i2 = indices[i], indices[i + 1], indices[i + 2]
		if type(i0) ~= "number" or type(i1) ~= "number" or type(i2) ~= "number" or
			i0 ~= math.floor(i0) or i1 ~= math.floor(i1) or i2 ~= math.floor(i2) or
			i0 < 0 or i0 >= vertex_count or
			i1 < 0 or i1 >= vertex_count or
			i2 < 0 or i2 >= vertex_count then
			mesh_x, mesh_y, mesh_tris = {}, {}, {}
			return false
		end
		mesh_tris[#mesh_tris + 1] = { i0 + 1, i1 + 1, i2 + 1 }
	end

	mesh_n_verts = vertex_count
	return true
end

local generated = false
local module_error = nil
if not module_ok then
	module_error = "obj.module(): " .. tostring(mesh_module)
elseif type(mesh_module) ~= "table" or type(mesh_module.generate) ~= "function" then
	module_error = "generate関数が登録されていません"
else
	local call_ok, vertices, indices = pcall(function()
		local data, pw, ph = obj.getpixeldata("object", "rgba")
		return mesh_module.generate(data, pw, ph, threshold, density, border)
	end)
	if not call_ok then
		module_error = "generate(): " .. tostring(vertices)
	else
		generated = import_module_mesh(vertices, indices)
		if not generated then module_error = "不正または空のメッシュが返されました" end
	end
end

if module_error and not puppet_geometry_error_reported then
	print("@warn", "puppet_geometry.mod2: " .. module_error .. " (格子メッシュを使用します)")
	puppet_geometry_error_reported = true
elseif generated then
	puppet_geometry_error_reported = false
end

if not generated then
	-- モジュールが無い、エラー、または空メッシュの場合の格子フォールバック。
	mesh_x, mesh_y, mesh_tris = {}, {}, {}
	local cell_size = math.max(w, h) / density
	local cols = math.max(1, math.ceil(w / cell_size))
	local rows = math.max(1, math.ceil(h / cell_size))
	local cw = w / cols
	local ch = h / rows

	mesh_n_verts = (rows + 1) * (cols + 1)
	for iy = 0, rows do
		for ix = 0, cols do
			local idx = iy * (cols + 1) + ix + 1
			mesh_x[idx] = ix * cw - hw
			mesh_y[idx] = iy * ch - hh
		end
	end

	local function grid_idx(ix, iy)
		return iy * (cols + 1) + ix + 1
	end

	for iy = 0, rows - 1 do
		for ix = 0, cols - 1 do
			local i0 = grid_idx(ix,     iy)
			local i1 = grid_idx(ix + 1, iy)
			local i2 = grid_idx(ix + 1, iy + 1)
			local i3 = grid_idx(ix,     iy + 1)

			if check_alpha_at((mesh_x[i0]+mesh_x[i1]+mesh_x[i2])/3, (mesh_y[i0]+mesh_y[i1]+mesh_y[i2])/3) or
				check_alpha_at(mesh_x[i0], mesh_y[i0]) or
				check_alpha_at(mesh_x[i1], mesh_y[i1]) or
				check_alpha_at(mesh_x[i2], mesh_y[i2]) then
				table.insert(mesh_tris, { i0, i1, i2 })
			end
			if check_alpha_at((mesh_x[i0]+mesh_x[i2]+mesh_x[i3])/3, (mesh_y[i0]+mesh_y[i2]+mesh_y[i3])/3) or
				check_alpha_at(mesh_x[i0], mesh_y[i0]) or
				check_alpha_at(mesh_x[i2], mesh_y[i2]) or
				check_alpha_at(mesh_x[i3], mesh_y[i3]) then
				table.insert(mesh_tris, { i0, i2, i3 })
			end
		end
	end
end
