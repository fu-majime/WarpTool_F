-- ベンド/詳細/スターチピンを、変形エンジンが扱える位置拘束へ展開する。
-- 変形前後のメッシュからピン位置の自然な局所姿勢を推定するため、
-- ベンドとスターチは位置ピンに追従しつつ周辺の形だけを制御できる。
local puppet_pin_dynamics = {}

local EPSILON = 1e-8

local function append_constraint(sources, destinations, layers, ranges,
		sx, sy, dx, dy, layer, range)
	local n = #sources
	sources[n + 1], sources[n + 2] = sx, sy
	destinations[n + 1], destinations[n + 2] = dx, dy
	layers[#layers + 1] = layer or 0
	ranges[#ranges + 1] = range or 0
end

-- 距離の近い頂点を使った重み付きsimilarity fitting。
-- 戻り値はピン中心の自然な移動先、回転、拡大率。
local function estimate_pose(vertices, deformed, px, py, radius)
	local samples = {}
	for i = 1, #vertices, 2 do
		local dx, dy = vertices[i] - px, vertices[i + 1] - py
		samples[#samples + 1] = { i = i, distance2 = dx * dx + dy * dy }
	end
	table.sort(samples, function(a, b) return a.distance2 < b.distance2 end)

	local count = math.min(12, #samples)
	if count == 0 then return px, py, 0, 1 end
	local floor_distance2 = math.max(radius * radius * 0.04, EPSILON)
	local weight_sum, spx, spy, dpx, dpy = 0, 0, 0, 0, 0
	for n = 1, count do
		local sample = samples[n]
		local weight = 1 / math.max(sample.distance2, floor_distance2)
		local i = sample.i
		weight_sum = weight_sum + weight
		spx, spy = spx + vertices[i] * weight, spy + vertices[i + 1] * weight
		dpx, dpy = dpx + deformed[i] * weight, dpy + deformed[i + 1] * weight
	end
	spx, spy = spx / weight_sum, spy / weight_sum
	dpx, dpy = dpx / weight_sum, dpy / weight_sum

	local real, imaginary, denominator = 0, 0, 0
	for n = 1, count do
		local sample = samples[n]
		local weight = 1 / math.max(sample.distance2, floor_distance2)
		local i = sample.i
		local sx, sy = vertices[i] - spx, vertices[i + 1] - spy
		local dx, dy = deformed[i] - dpx, deformed[i + 1] - dpy
		real = real + weight * (sx * dx + sy * dy)
		imaginary = imaginary + weight * (sx * dy - sy * dx)
		denominator = denominator + weight * (sx * sx + sy * sy)
	end

	local angle, scale = 0, 1
	if denominator > EPSILON then
		angle = math.atan2(imaginary, real)
		scale = math.sqrt(real * real + imaginary * imaginary) / denominator
		-- 壊れかけたメッシュや画像外ピンで極端な拘束を作らない。
		scale = math.max(0.2, math.min(5, scale))
	end
	local cosine, sine = math.cos(angle) * scale, math.sin(angle) * scale
	local ox, oy = px - spx, py - spy
	local target_x = dpx + ox * cosine - oy * sine
	local target_y = dpy + ox * sine + oy * cosine
	return target_x, target_y, angle, scale
end

-- 各頂点を最寄りのピンに割り当て、その領域内から8方向の代表点を拾う。
-- 半径は隣のピンまでの距離を基準にするが、メッシュ全体に対して上限を
-- 設ける。最遠点を直接使うと、疎なピン配置で腕などを巻き込んで長大な
-- 三角形を作るためである。
local function handle_footprint(vertices, pin_data, pin, minimum_radius)
	local footprint = {}
	for i = 1, #vertices, 2 do
		local ox, oy = vertices[i] - pin.sx, vertices[i + 1] - pin.sy
		local own_distance2 = ox * ox + oy * oy
		local owned = true
		for _, other in ipairs(pin_data) do
			if other ~= pin then
				local dx, dy = vertices[i] - other.sx, vertices[i + 1] - other.sy
				if dx * dx + dy * dy + EPSILON < own_distance2 then
					owned = false
					break
				end
			end
		end
		if owned and own_distance2 > EPSILON then
			footprint[#footprint + 1] = { x = vertices[i], y = vertices[i + 1] }
		end
	end
	if #footprint == 0 then
		for sector = 0, 7 do
			local angle = sector / 8 * math.pi * 2
			footprint[#footprint + 1] = {
				x = pin.sx + math.cos(angle) * minimum_radius,
				y = pin.sy + math.sin(angle) * minimum_radius,
			}
		end
	end
	return footprint
end

local function append_handle(sources, destinations, layers, ranges, pin_data, vertices, pin,
		tx, ty, angle, scale, minimum_radius)
	local px, py, layer = pin.sx, pin.sy, pin.layer
	append_constraint(sources, destinations, layers, ranges, px, py, tx, ty, layer, 0)
	local cosine, sine = math.cos(angle) * scale, math.sin(angle) * scale
	local footprint = handle_footprint(vertices, pin_data, pin, minimum_radius)
	local debug_points = {}
	for _, point in ipairs(footprint) do
		local ox, oy = point.x - px, point.y - py
		local dx = tx + ox * cosine - oy * sine
		local dy = ty + ox * sine + oy * cosine
		append_constraint(sources, destinations, layers, ranges,
			point.x, point.y, dx, dy, layer, 0)
		debug_points[#debug_points + 1] = {
			sx = point.x, sy = point.y, dx = dx, dy = dy,
		}
	end
	return debug_points
end

function puppet_pin_dynamics.build(vertices, preliminary, pin_data, radius, pin_type)
	local sources, destinations, layers, ranges, debug_handles = {}, {}, {}, {}, {}
	for _, pin in ipairs(pin_data) do
		local kind = pin.kind
		if kind == pin_type.POSITION or kind == pin_type.BONE then
			append_constraint(sources, destinations, layers, ranges,
				pin.sx, pin.sy, pin.dx, pin.dy, pin.layer, 0)
		elseif kind == pin_type.DETAIL then
			local _, _, natural_angle, natural_scale = estimate_pose(
				vertices, preliminary, pin.sx, pin.sy, radius)
			local points = append_handle(sources, destinations, layers, ranges,
				pin_data, vertices, pin,
				pin.dx, pin.dy,
				natural_angle + math.rad(pin.rotation), natural_scale * pin.scale,
				radius)
			debug_handles[#debug_handles + 1] = {
				kind = kind, sx = pin.sx, sy = pin.sy, dx = pin.dx, dy = pin.dy,
				points = points,
			}
		elseif kind == pin_type.BEND then
			local tx, ty, natural_angle, natural_scale = estimate_pose(
				vertices, preliminary, pin.sx, pin.sy, radius)
			local points = append_handle(sources, destinations, layers, ranges,
				pin_data, vertices, pin,
				tx, ty,
				natural_angle + math.rad(pin.rotation), natural_scale * pin.scale,
				radius)
			debug_handles[#debug_handles + 1] = {
				kind = kind, sx = pin.sx, sy = pin.sy, dx = tx, dy = ty,
				points = points,
			}
		elseif kind == pin_type.STARCH then
			local tx, ty, natural_angle = estimate_pose(
				vertices, preliminary, pin.sx, pin.sy, radius)
			-- 自然な移動と回転には追従するが、局所拡縮は捨てて形を保つ。
			local points = append_handle(sources, destinations, layers, ranges,
				pin_data, vertices, pin,
				tx, ty, natural_angle, 1, radius)
			debug_handles[#debug_handles + 1] = {
				kind = kind, sx = pin.sx, sy = pin.sy, dx = tx, dy = ty,
				points = points,
			}
		elseif kind == pin_type.OVERLAP then
			local tx, ty = estimate_pose(vertices, preliminary, pin.sx, pin.sy, radius)
			-- 負値はRust側で「変形拘束ではない重なり専用ピン」と識別する。
			-- -1なら実範囲0、以降は -(range + 1) で衝突しない。
			append_constraint(sources, destinations, layers, ranges,
				pin.sx, pin.sy, tx, ty, pin.layer, -(pin.range + 1))
			debug_handles[#debug_handles + 1] = {
				kind = kind, sx = pin.sx, sy = pin.sy, dx = tx, dy = ty,
				points = {},
			}
		end
	end
	return sources, destinations, layers, ranges, debug_handles
end
