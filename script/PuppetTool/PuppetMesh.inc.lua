-- メッシュ生成
local mesh_x, mesh_y = {}, {}
local mesh_tris = {}
local mesh_n_verts = 0

local function import_module_mesh(vertices, indices)
	if type(vertices) == "table" and type(indices) == "table" and #vertices == 0 and #indices == 0 then return true end
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

if module_error then
	error("puppet_geometry.mod2: " .. module_error)
end
