-- 方向波形ワープ共通UI
--group:波形,true
---$select:波形
---サイン=1
---矩形=2
---三角形=3
---のこぎり(正)=4
---のこぎり(負)=5
---円形=6
---半円形(正)=7
---半円形(負)=8
local shape = 1

---$check:鏡面反射
local mirror = 0

--group:ランダム,false
---$track:高さランダム[%], min = 0, max = 100
local randHeight = 0
---$track:幅ランダム[%], min = 0, max = 100
local randWidth = 0
---$track:ランダムシード, min = -65536, max = 65535, step = 1
local randSeed = 0

--group:その他,false
---$select:固定処理
---なし=1
---すべて=2
---上下=3
---左右=4
---上=5
---下=6
---左=7
---右=8
local fix = 1

---$track:フェード幅[px], min = 0, max = 1000
local fadePixel = 20
---$value:PI
local PI = {}
