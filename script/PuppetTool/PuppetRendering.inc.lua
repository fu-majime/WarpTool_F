-- =============================================================
-- 5. 描画バッファ設定
-- =============================================================
local dstw, dsth = w, h
local sx_off, sy_off = 0, 0

if obj.getinfo("filter") then
	obj.setoption("drawtarget", "tempbuffer", dstw, dsth)
else
	local minbx, minby = 1e9, 1e9
	local maxbx, maxby = -1e9, -1e9
	for i = 1, mesh_n_verts do
		if def_x[i] and def_x[i] < minbx then minbx = def_x[i] end
		if def_x[i] and def_x[i] > maxbx then maxbx = def_x[i] end
		if def_y[i] and def_y[i] < minby then minby = def_y[i] end
		if def_y[i] and def_y[i] > maxby then maxby = def_y[i] end
	end
	-- ガイドも描画バッファへ収め、長く移動したピンの線や点を欠けさせない。
	if show_gui then
		for i = 1, pins do
			minbx = math.min(minbx, pin_sx[i], pin_tx[i])
			maxbx = math.max(maxbx, pin_sx[i], pin_tx[i])
			minby = math.min(minby, pin_sy[i], pin_ty[i])
			maxby = math.max(maxby, pin_sy[i], pin_ty[i])
		end
		for _, handle in ipairs(pin_dynamics_debug) do
			minbx, maxbx = math.min(minbx, handle.dx), math.max(maxbx, handle.dx)
			minby, maxby = math.min(minby, handle.dy), math.max(maxby, handle.dy)
			for _, point in ipairs(handle.points) do
				minbx = math.min(minbx, point.sx, point.dx)
				maxbx = math.max(maxbx, point.sx, point.dx)
				minby = math.min(minby, point.sy, point.dy)
				maxby = math.max(maxby, point.sy, point.dy)
			end
		end
	end
	if minbx > maxbx then minbx = -hw; maxbx = hw end
	if minby > maxby then minby = -hh; maxby = hh end
	local mx = (maxbx - minbx) * 0.1
	local my = (maxby - minby) * 0.1
	if show_gui then
		mx = math.max(mx, 16)
		my = math.max(my, 16)
	end
	minbx = minbx - mx; maxbx = maxbx + mx
	minby = minby - my; maxby = maxby + my
	dstw = math.max(1, math.ceil(maxbx - minbx))
	dsth = math.max(1, math.ceil(maxby - minby))
	obj.setoption("drawtarget", "tempbuffer", dstw, dsth)
	sx_off = (minbx + maxbx) / 2
	sy_off = (minby + maxby) / 2
end

-- =============================================================
-- 6. Rust側で細分化・Zソート済みの三角形を描画
-- =============================================================
local flat_vtx = {}
if #render_vertices % 12 ~= 0 then
	error("puppet_geometry: 描画頂点配列の長さが不正です (" .. #render_vertices .. ")")
end
for i = 1, #render_vertices, 4 do
	flat_vtx[#flat_vtx + 1] = {
		render_vertices[i] - sx_off, render_vertices[i + 1] - sy_off, 0,
		render_vertices[i + 2], render_vertices[i + 3]
	}
end
if #flat_vtx >= 3 then
	obj.setoption("blend", "alpha_add")
	obj.drawpoly(flat_vtx, 3)
end

-- =============================================================
-- 7. 可視化（show=true）
-- =============================================================
if show_gui then
	local COLORS = {
		pink   = 0xff69b4,
		orange = 0xffa020,
		red    = 0xe02848,
		green  = 0x20e060,
		blue   = 0x4080ff,
		purple = 0xa060e0,
		lime   = 0x80ff40,
		gray   = 0x808080,
		cyan   = 0x00ccff,
		yellow = 0xf2bf26
	}

	local COLOR_THEME = {
		mesh           = { color = COLORS.cyan,   alpha = 0.45 },
		arrowPath      = { color = COLORS.gray,   alpha = 0.80 },
		arrowProgress  = { color = COLORS.lime,   alpha = 0.90 },
		pinSource      = {
			[PIN_TYPE.POSITION] = { color = COLORS.yellow, alpha = 0.55 },
			[PIN_TYPE.BONE]     = { color = COLORS.green,  alpha = 0.55 },
			[PIN_TYPE.BEND]     = { color = COLORS.orange, alpha = 0.55 },
			[PIN_TYPE.DETAIL]   = { color = COLORS.red,    alpha = 0.55 },
			[PIN_TYPE.OVERLAP]  = { color = COLORS.blue,   alpha = 0.55 },
			[PIN_TYPE.STARCH]   = { color = COLORS.purple, alpha = 0.55 }
		},
		pinDestination = { color = COLORS.pink, alpha = 0.75 }
	}

	local function premultiplied_rgba(theme)
		local rgb = theme.color
		local r = math.floor(rgb / 0x10000) / 0xff
		local g = math.floor(rgb / 0x100) % 0x100 / 0xff
		local b = rgb % 0x100 / 0xff
		local a = theme.alpha
		return r * a, g * a, b * a, a
	end

	local function overlap_visual(x, y)
		local influence = 0
		for _, pin in ipairs(pin_data) do
			if pin.kind == PIN_TYPE.OVERLAP and pin.show_range
				and pin.range > 0 and pin.layer ~= 0 then
				local dx, dy = x - pin.sx, y - pin.sy
				local distance = math.sqrt(dx * dx + dy * dy)
				if distance < pin.range then
					local t = math.max(0, math.min(1, 1 - distance / pin.range))
					local falloff = t * t * (3 - 2 * t)
					influence = influence + pin.layer * falloff
				end
			end
		end
		return math.max(-1, math.min(1, influence))
	end

	-- 変形後メッシュと同じ三角形に白黒半透明面を貼る。
	-- ワイヤー自体の色は変えず、重なりの強さを面のアルファで示す。
	local overlap_plane = {}
	for _, tri in ipairs(mesh_tris) do
		local triangle_vertices, visible = {}, false
		for vertex = 1, 3 do
			local index = tri[vertex]
			local influence = overlap_visual(mesh_x[index], mesh_y[index])
			local alpha = math.abs(influence) * 0.65
			visible = visible or alpha > 0.001
			local color = influence > 0 and alpha or 0
			triangle_vertices[vertex] = {
				def_x[index] - sx_off, def_y[index] - sy_off, 0,
				color, color, color, alpha,
			}
		end
		if visible then
			for vertex = 1, 3 do
				overlap_plane[#overlap_plane + 1] = triangle_vertices[vertex]
			end
		end
	end
	if #overlap_plane >= 3 then
		obj.setoption("blend", "alpha_add")
		obj.drawpoly(overlap_plane, 3)
	end

	-- ワイヤーフレーム
	local wire = {}
	local lw = 1.2
	local cr, cg, cb, ca = premultiplied_rgba(COLOR_THEME.mesh)
	local drawn_edges = {}

	for ti = 1, #mesh_tris do
		local tri = mesh_tris[ti]
		local vx = {
			def_x[tri[1]] - sx_off, def_y[tri[1]] - sy_off,
			def_x[tri[2]] - sx_off, def_y[tri[2]] - sy_off,
			def_x[tri[3]] - sx_off, def_y[tri[3]] - sy_off
		}
		for e = 0, 2 do
			local ne = ((e + 1) % 3)
			local ia, ib = tri[e + 1], tri[ne + 1]
			local edge_key = math.min(ia, ib) .. ":" .. math.max(ia, ib)
			if not drawn_edges[edge_key] then
				drawn_edges[edge_key] = true
				local ax, ay = vx[e * 2 + 1], vx[e * 2 + 2]
				local bx, by = vx[ne * 2 + 1], vx[ne * 2 + 2]
				local edx, edy = bx - ax, by - ay
				local elen = math.sqrt(edx * edx + edy * edy)
				if elen > 0.1 then
					local nx = -edy / elen * lw
					local ny =  edx / elen * lw
					table.insert(wire, { ax + nx, ay + ny, 0, cr, cg, cb, ca })
					table.insert(wire, { bx + nx, by + ny, 0, cr, cg, cb, ca })
					table.insert(wire, { bx - nx, by - ny, 0, cr, cg, cb, ca })
					table.insert(wire, { ax - nx, ay - ny, 0, cr, cg, cb, ca })
				end
			end
		end
	end

	if #wire >= 4 then
		obj.setoption("blend", "alpha_add")
		obj.drawpoly(wire)
	end

	-- ボーンは親側を幅広く、子側を尖らせた平面三角形で表示する。
	local bone_faces = {}
	local bone_r, bone_g, bone_b, bone_a = premultiplied_rgba {
		color = COLORS.green, alpha = 0.42
	}
	local bone_forest_for_render = bone_forest
	while type(bone_forest_for_render) == "table" and #bone_forest_for_render == 1
		and type(bone_forest_for_render[1]) == "table"
		and type(bone_forest_for_render[1][1]) == "table" do
		bone_forest_for_render = bone_forest_for_render[1]
	end
	local function append_bone(tree, parent_id)
		if type(tree) ~= "table" then return end
		local pin_id = math.floor(tonumber(tree[1]) or 0)
		if parent_id and pin_id >= 1 and pin_id <= pins then
			local px, py = pin_dx[parent_id] - sx_off, pin_dy[parent_id] - sy_off
			local cx, cy = pin_dx[pin_id] - sx_off, pin_dy[pin_id] - sy_off
			local dx, dy = cx - px, cy - py
			local length = math.sqrt(dx * dx + dy * dy)
			if length > 0.1 then
				local half_width = math.min(14, length * 0.28)
				local nx, ny = -dy / length * half_width, dx / length * half_width
				bone_faces[#bone_faces + 1] = { px + nx, py + ny, 0, bone_r, bone_g, bone_b, bone_a }
				bone_faces[#bone_faces + 1] = { px - nx, py - ny, 0, bone_r, bone_g, bone_b, bone_a }
				bone_faces[#bone_faces + 1] = { cx, cy, 0, bone_r, bone_g, bone_b, bone_a }
			end
		end
		for i = 2, #tree do append_bone(tree[i], pin_id) end
	end
	if type(bone_forest_for_render) == "table" then
		for i = 1, #bone_forest_for_render do append_bone(bone_forest_for_render[i], nil) end
	end
	if #bone_faces >= 3 then
		obj.setoption("blend", "alpha_add")
		obj.drawpoly(bone_faces, 3)
	end

	local function make_line(x1, y1, x2, y2, width, theme)
		local dx, dy = x2 - x1, y2 - y1
		local len = math.sqrt(dx * dx + dy * dy)
		if len <= 0.1 then return nil end
		local nx, ny = -dy / len * width * 0.5, dx / len * width * 0.5
		local r, g, b, a = premultiplied_rgba(theme)
		return {
			x1 + nx, y1 + ny, 0, x2 + nx, y2 + ny, 0,
			x2 - nx, y2 - ny, 0, x1 - nx, y1 - ny, 0,
			r, g, b, a, r, g, b, a, r, g, b, a, r, g, b, a
		}
	end

	-- ベンド/詳細/スターチの自動追従を可視化する。
	-- 中心線がピンの追従移動、薄い放射線が外周拘束の移動、
	-- 色付きの輪郭が実際に回転・拡縮へ使った領域を表す。
	local dynamics_lines = {}
	for _, handle in ipairs(pin_dynamics_debug) do
		if handle.kind ~= PIN_TYPE.STARCH and handle.kind ~= PIN_TYPE.OVERLAP then
			local source_theme = COLOR_THEME.pinSource[handle.kind] or COLOR_THEME.arrowProgress
			local influence_theme = { color = source_theme.color, alpha = 0.38 }
			local center_line = make_line(
				handle.sx - sx_off, handle.sy - sy_off,
				handle.dx - sx_off, handle.dy - sy_off, 4, source_theme)
			if center_line then dynamics_lines[#dynamics_lines + 1] = center_line end
			for index, point in ipairs(handle.points) do
				local movement = make_line(
					point.sx - sx_off, point.sy - sy_off,
					point.dx - sx_off, point.dy - sy_off, 1.5, influence_theme)
				if movement then dynamics_lines[#dynamics_lines + 1] = movement end
				if #handle.points > 1 then
					local next_point = handle.points[index % #handle.points + 1]
					local boundary = make_line(
						point.dx - sx_off, point.dy - sy_off,
						next_point.dx - sx_off, next_point.dy - sy_off, 2, influence_theme)
					if boundary then dynamics_lines[#dynamics_lines + 1] = boundary end
				end
			end
		end
	end
	if #dynamics_lines > 0 then
		obj.setoption("blend", "alpha_add")
		obj.drawpoly(dynamics_lines)
	end

	-- src→dstの全経路と、src→現在位置（ratio反映）の矢印軸を一括描画する。
	local pin_lines = {}
	local arrow_lines = {}
	local arrow_heads = {}
	for i = 1, pins do
		local kind = pin_types[i]
		if kind ~= PIN_TYPE.BONE and kind ~= PIN_TYPE.STARCH
			and kind ~= PIN_TYPE.OVERLAP then
		local ox, oy = pin_sx[i] - sx_off, pin_sy[i] - sy_off
		local fx, fy = pin_tx[i] - sx_off, pin_ty[i] - sy_off
		local cx, cy = pin_dx[i] - sx_off, pin_dy[i] - sy_off
		local line = make_line(ox, oy, fx, fy, 8, COLOR_THEME.arrowPath)
		if line then table.insert(pin_lines, line) end

		local adx, ady = cx - ox, cy - oy
		local alen = math.sqrt(adx * adx + ady * ady)
		if alen > 1 then
			local arrow = make_line(ox, oy, cx, cy, 8, COLOR_THEME.arrowProgress)
			if arrow then table.insert(arrow_lines, arrow) end

			-- 三角形図形の上端中央が現在位置に来るよう、UV付き四角形を置く。
			local hx, hy = adx / alen, ady / alen
			local nx, ny = -hy, hx
			local head_len = 24 * math.max(math.min(1, (ratio/100)/0.1, (1-ratio/100)/0.1), 0)
			local head_width = head_len * 1.5
			local bx, by = cx - hx * head_len, cy - hy * head_len
			table.insert(arrow_heads, { cx - nx * head_width * 0.5, cy - ny * head_width * 0.5, 0, 0, 0 })
			table.insert(arrow_heads, { cx + nx * head_width * 0.5, cy + ny * head_width * 0.5, 0, 1, 0 })
			table.insert(arrow_heads, { bx + nx * head_width * 0.5, by + ny * head_width * 0.5, 0, 1, 1 })
			 table.insert(arrow_heads, { bx - nx * head_width * 0.5, by - ny * head_width * 0.5, 0, 0, 1 })
		end
		end
	end
	if #pin_lines > 0 then
		obj.setoption("blend", "alpha_add")
		obj.drawpoly(pin_lines)
	end
	if #arrow_lines > 0 then
		obj.drawpoly(arrow_lines)
	end

	-- 矢じりと点は既成図形を使う。操作側の点だけ大きく、不透明にする。
	local saved_props = {
		obj.ox, obj.oy, obj.oz, obj.cx, obj.cy, obj.cz,
		obj.rx, obj.ry, obj.rz, obj.sx, obj.sy, obj.sz, obj.alpha
	}
	obj.setoption("drawtarget", "tempbuffer")
	local point_base_size = 32
	local function draw_destination_points(xs, ys, theme, active)
		if not obj.load("figure", "円", theme.color, point_base_size) then return end
		local size = active and 32 or 16
		local alpha = active and 1.0 or theme.alpha
		for i = 1, pins do
			local kind = pin_types[i]
			if kind ~= PIN_TYPE.BONE and kind ~= PIN_TYPE.BEND and kind ~= PIN_TYPE.STARCH
				and kind ~= PIN_TYPE.OVERLAP then
				obj.draw(xs[i] - sx_off, ys[i] - sy_off, 0, size / point_base_size, alpha)
			end
		end
	end
	local function draw_source_points(active)
		local size = active and 32 or 16
		for pin_type = PIN_TYPE.POSITION, PIN_TYPE.OVERLAP do
			local theme = COLOR_THEME.pinSource[pin_type]
			if theme and obj.load("figure", "円", theme.color, point_base_size) then
				local alpha = active and 1.0 or theme.alpha
				for i = 1, pins do
					if pin_types[i] == pin_type and pin_type ~= PIN_TYPE.BONE
						and pin_type ~= PIN_TYPE.STARCH and pin_type ~= PIN_TYPE.OVERLAP then
						obj.draw(pin_sx[i] - sx_off, pin_sy[i] - sy_off, 0,
							size / point_base_size, alpha)
					end
				end
			end
		end
	end
	draw_source_points(edit_source)
	draw_destination_points(pin_tx, pin_ty, COLOR_THEME.pinDestination, not edit_source)
	-- ボーンピンは移動元を描かず、現在位置を菱形で示す。
	local diamond_vertices = {}
	local diamond_size = edit_source and 9 or 16
	local dr, dg, db, da = premultiplied_rgba {
		color = COLORS.green, alpha = edit_source and 0.7 or 1.0
	}
	for i = 1, pins do
		if pin_types[i] == PIN_TYPE.BONE then
			local x, y = pin_dx[i] - sx_off, pin_dy[i] - sy_off
			diamond_vertices[#diamond_vertices + 1] = { x, y - diamond_size, 0, dr, dg, db, da }
			diamond_vertices[#diamond_vertices + 1] = { x + diamond_size, y, 0, dr, dg, db, da }
			diamond_vertices[#diamond_vertices + 1] = { x, y + diamond_size, 0, dr, dg, db, da }
			diamond_vertices[#diamond_vertices + 1] = { x - diamond_size, y, 0, dr, dg, db, da }
		end
	end
	if #diamond_vertices >= 4 then obj.drawpoly(diamond_vertices) end
	-- 自動追従した中心は、ピン種別色の小さい不透明点で示す。
	for pin_type = PIN_TYPE.BEND, PIN_TYPE.OVERLAP do
		local theme = COLOR_THEME.pinSource[pin_type]
		if theme and obj.load("figure", "円", theme.color, point_base_size) then
			for _, handle in ipairs(pin_dynamics_debug) do
				if handle.kind == pin_type then
					obj.draw(handle.dx - sx_off, handle.dy - sy_off, 0,
						12 / point_base_size, 1)
				end
			end
		end
	end
	if #arrow_heads >= 4 and obj.load("figure", "三角形", COLOR_THEME.arrowProgress.color, 16) then
		obj.drawpoly(arrow_heads, 4, COLOR_THEME.arrowProgress.alpha)
	end
	obj.ox, obj.oy, obj.oz, obj.cx, obj.cy, obj.cz,
	obj.rx, obj.ry, obj.rz, obj.sx, obj.sy, obj.sz, obj.alpha = unpack(saved_props)
end
