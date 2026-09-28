local puppet_pin_dynamics = {}
local EPSILON = 1e-8

function puppet_pin_dynamics.bend_frames(pin_data, pin_type)
    local found=false
    for _,pin in ipairs(pin_data) do
        if pin.kind==pin_type.BEND then found=true
        elseif pin.kind~=pin_type.OVERLAP then return false end
    end
    return found
end

local function append(sources, destinations, layers, ranges, sx, sy, dx, dy, layer, range)
    sources[#sources+1], sources[#sources+2] = sx, sy
    destinations[#destinations+1], destinations[#destinations+2] = dx, dy
    layers[#layers+1], ranges[#ranges+1] = layer or 0, range or 0
end

local function topology(vertices, indices)
    local graph = {}
    for i=1,#vertices/2 do graph[i] = {} end
    for i=1,#indices,3 do
        for k=0,2 do
            local a,b = indices[i+k]+1,indices[i+(k+1)%3]+1
            graph[a][b],graph[b][a] = true,true
        end
    end
    return graph
end

local function neighborhood(vertices, graph, pin, radius)
    local center, best = 1,math.huge
    for i=1,#vertices/2 do
        local d = (vertices[i*2-1]-pin.sx)^2+(vertices[i*2]-pin.sy)^2
        if d<best then best,center = d,i end
    end
    local distances, queue = {[center]=0},{center}
    local cursor = 1
    while cursor<=#queue do
        local a=queue[cursor]; cursor=cursor+1
        for b in pairs(graph[a]) do
            local d=distances[a]+math.sqrt((vertices[a*2-1]-vertices[b*2-1])^2+(vertices[a*2]-vertices[b*2])^2)
            if d<radius*2 and (not distances[b] or d<distances[b]-EPSILON) then
                distances[b]=d; queue[#queue+1]=b
            end
        end
    end
    local samples={}
    for i,d in pairs(distances) do samples[#samples+1]={vertex=i,distance=d} end
    table.sort(samples,function(a,b)
        if a.distance==b.distance then return a.vertex<b.vertex end
        return a.distance<b.distance
    end)
    return samples,center,best
end

function puppet_pin_dynamics.overlap_influence(vertices, indices, pin_data, overlap_type)
    local graph=topology(vertices,indices)
    local influence={}
    for i=1,#vertices/2 do influence[i]=0 end
    for _,pin in ipairs(pin_data) do
        if pin.kind==overlap_type and pin.show_range and (pin.range or 0)>0
            and (pin.layer or 0)~=0 then
            local samples=neighborhood(vertices,graph,pin,pin.range)
            for _,sample in ipairs(samples) do
                if sample.distance<pin.range then
                    local t=math.max(0,math.min(1,1-sample.distance/pin.range))
                    local falloff=t*t*(3-2*t)
                    influence[sample.vertex]=influence[sample.vertex]+pin.layer*falloff
                end
            end
        end
    end
    for i=1,#influence do influence[i]=math.max(-1,math.min(1,influence[i])) end
    return influence
end

local function natural_pose(vertices, deformed, pin, samples, center, radius)
    local real, imaginary, denominator = 0,0,0
    local px,py = vertices[center*2-1],vertices[center*2]
    local qx,qy = deformed[center*2-1],deformed[center*2]
    -- Fit the derivative about the actual material point, not a drifting
    -- centroid of twelve arbitrarily close vertices on another limb.
    for _,sample in ipairs(samples) do
        local i=sample.vertex*2-1
        local x,y=vertices[i]-px,vertices[i+1]-py
        local u,v=deformed[i]-qx,deformed[i+1]-qy
        local weight=1/math.max(sample.distance^2,radius^2*0.25)
        real=real+weight*(x*u+y*v)
        imaginary=imaginary+weight*(x*v-y*u)
        denominator=denominator+weight*(x*x+y*y)
    end
    local angle,scale=0,1
    if denominator>EPSILON then
        angle=math.atan2(imaginary,real)
        scale=math.sqrt(real*real+imaginary*imaginary)/denominator
    end
    local c,s=math.cos(angle)*scale,math.sin(angle)*scale
    return {x=qx+(pin.sx-px)*c-(pin.sy-py)*s,
        y=qy+(pin.sx-px)*s+(pin.sy-py)*c,rotation=angle,scale=scale}
end

local function material_radius(vertices, indices, pin, fallback)
    local edges={}
    for i=1,#indices,3 do
        for k=0,2 do
            local a,b=indices[i+k]+1,indices[i+(k+1)%3]+1
            local key=math.min(a,b)..":"..math.max(a,b)
            local e=edges[key] or {a=a,b=b,n=0}; e.n=e.n+1; edges[key]=e
        end
    end
    local best=math.huge
    for _,e in pairs(edges) do
        if e.n==1 then
            local ax,ay=vertices[e.a*2-1],vertices[e.a*2]
            local x,y=vertices[e.b*2-1]-ax,vertices[e.b*2]-ay
            local t=math.max(0,math.min(1,((pin.sx-ax)*x+(pin.sy-ay)*y)/math.max(x*x+y*y,EPSILON)))
            best=math.min(best,math.sqrt((pin.sx-ax-t*x)^2+(pin.sy-ay-t*y)^2))
        end
    end
    return math.max(fallback,math.min(best,fallback*8))
end

function puppet_pin_dynamics.build(vertices, preliminary, pin_data, radius, pin_type, indices)
    local graph=topology(vertices,indices)
    local bend_frames=puppet_pin_dynamics.bend_frames(pin_data,pin_type)
    local sources,destinations,layers,ranges,debug,payload={},{},{},{},{},{}
    for _,pin in ipairs(pin_data) do
        local r=material_radius(vertices,indices,pin,radius)
        local samples,center=neighborhood(vertices,graph,pin,r)
        local program=natural_pose(vertices,preliminary,pin,samples,center,r)
        local moving=pin.kind==pin_type.POSITION or pin.kind==pin_type.DETAIL or pin.kind==pin_type.BONE
        local controlled=pin.kind==pin_type.DETAIL or pin.kind==pin_type.BEND
        local rotation=controlled and math.rad(pin.rotation or 0) or 0
        local scale=controlled and (pin.scale or 1) or 1
        local active=moving or pin.kind==pin_type.STARCH or (bend_frames and pin.kind==pin_type.BEND) or math.abs(rotation)>EPSILON or math.abs(scale-1)>EPSILON
        local explicit_rotation=pin.kind==pin_type.DETAIL or (pin.kind==pin_type.BEND and (bend_frames or math.abs(rotation)>EPSILON))
        local x,y=moving and pin.dx or program.x,moving and pin.dy or program.y
        if moving then append(sources,destinations,layers,ranges,pin.sx,pin.sy,x,y,pin.layer,0)
        elseif pin.kind==pin_type.OVERLAP then
            append(sources,destinations,layers,ranges,pin.sx,pin.sy,x,y,pin.layer,-((pin.range or 0)+1))
        end
        if active then
            -- Scale is relative to REST geometry, never the distorted preview.
            -- Bone hierarchy controls destination positions only. At identical
            -- positions a bone and a position pin have identical material rules.
            for _,v in ipairs({pin.sx,pin.sy,x,y,explicit_rotation and 0 or program.rotation,1,
                rotation,scale,r,moving and 1 or 0,explicit_rotation and 1 or 0,1}) do
                payload[#payload+1]=v
            end
        end
        debug[#debug+1]={kind=pin.kind,sx=pin.sx,sy=pin.sy,dx=x,dy=y,points={},
            program=program,user={position=moving,rotation=rotation,scale=scale}}
    end
    return sources,destinations,layers,ranges,debug,payload
end
