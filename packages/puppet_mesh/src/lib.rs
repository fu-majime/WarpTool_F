/// puppet_mesh — 輪郭ベースメッシュ生成 DLL
///
/// AviUtl2 LuaJIT FFI から呼び出されるC互換FFIエクスポート関数を提供する。
/// RGBA画像データからα境界に基づいた高品質メッシュを生成する。

mod contour;
mod delaunay;
mod dilate;
mod poisson;
mod simplify;

/// FFIで返却されるメッシュ結果。
/// LuaJIT の ffi.cdef で同じレイアウトを定義する。
#[repr(C)]
pub struct MeshResult {
    /// 頂点座標配列 [x0, y0, x1, y1, ...] (画像中心原点, f32)
    pub vertices: *mut f32,
    /// 頂点数 (座標配列の要素数は vertex_count * 2)
    pub vertex_count: u32,
    /// 三角形インデックス配列 [i0, i1, i2, ...] (0-based, u32)
    pub indices: *mut u32,
    /// インデックス数 (三角形数 = index_count / 3)
    pub index_count: u32,
}

/// メッシュを生成する。
///
/// # パラメータ
/// - `data`: RGBA画像データへのポインタ (各ピクセル4バイト: R,G,B,A)
/// - `w`, `h`: 画像サイズ（ピクセル）
/// - `threshold`: αしきい値 (0-255)
/// - `density`: メッシュ密度 (5-50, 大きいほど細かい)
/// - `border_px`: 縁取り膨張ピクセル数 (0=膨張なし)
///
/// # 戻り値
/// ヒープ上に確保された MeshResult へのポインタ。
/// 使用後は `puppet_mesh_free()` で解放すること。
#[no_mangle]
pub extern "C" fn puppet_mesh_generate(
    data: *const u8,
    w: i32,
    h: i32,
    threshold: i32,
    density: i32,
    border_px: i32,
) -> *mut MeshResult {
    let result = std::panic::catch_unwind(|| {
        generate_impl(data, w, h, threshold, density, border_px)
    });

    match result {
        Ok(r) => r,
        Err(_) => {
            // パニック時は空の結果を返す
            Box::into_raw(Box::new(MeshResult {
                vertices: std::ptr::null_mut(),
                vertex_count: 0,
                indices: std::ptr::null_mut(),
                index_count: 0,
            }))
        }
    }
}

fn generate_impl(
    data: *const u8,
    w: i32,
    h: i32,
    threshold: i32,
    density: i32,
    border_px: i32,
) -> *mut MeshResult {
    if data.is_null() || w <= 0 || h <= 0 {
        return Box::into_raw(Box::new(MeshResult {
            vertices: std::ptr::null_mut(),
            vertex_count: 0,
            indices: std::ptr::null_mut(),
            index_count: 0,
        }));
    }

    let w = w as usize;
    let h = h as usize;
    let threshold = threshold.clamp(0, 255) as u8;
    let density = density.clamp(5, 200) as usize;
    let border = border_px.max(0) as f32;

    // メッシュ間隔と誤差許容値の計算
    let max_dim = w.max(h) as f32;
    let min_spacing = max_dim / density as f32;
    // Douglas-Peuckerの許容誤差
    let epsilon = min_spacing * 0.4;

    // DPで内側に食い込む最大量が epsilon なので、
    // 事前に border + epsilon だけ膨張させることで、
    // 単純化後のメッシュ辺が不透明領域を絶対に切り取らないようにする。
    let total_dilate = border + epsilon;

    // 画像の境界に接するオブジェクトも正しく輪郭 추출・膨張できるようにパディングする
    let pad = (total_dilate.ceil() as usize) + 2;
    let padded_w = w + 2 * pad;
    let padded_h = h + 2 * pad;

    // RGBA画像からα値を抽出し、パディングされた二値マップを作成
    let rgba = unsafe { std::slice::from_raw_parts(data, w * h * 4) };
    let mut binary = vec![false; padded_w * padded_h];
    for y in 0..h {
        let src_row_start = y * w * 4;
        let dst_row_start = (y + pad) * padded_w + pad;
        for x in 0..w {
            let a = rgba[src_row_start + x * 4 + 3];
            binary[dst_row_start + x] = a >= threshold;
        }
    }

    // EDT膨張
    if total_dilate > 0.0 {
        dilate::dilate_edt(&mut binary, padded_w, padded_h, total_dilate);
    }

    // 膨張後のパディング済みα画像（Delaunayの三角形カリング・Poisson用）
    let dilated_alpha: Vec<u8> = binary.iter().map(|&b| if b { 255 } else { 0 }).collect();

    // 1. 輪郭抽出 (Marching Squares, パディングされたマップから)
    // 座標はパディング済みグリッドの中心が原点になるが、左右上下対称にパディング
    // しているため、元の画像の中心と物理的に全く同じ位置になる！そのまま処理可能。
    let contour_loops = contour::extract_contours_from_binary(&binary, padded_w, padded_h);

    // 2. 輪郭単純化 (Douglas-Peucker)
    let mut simplified_contours: Vec<Vec<(f32, f32)>> = Vec::new();
    for loop_pts in &contour_loops {
        let simplified = if loop_pts.len() > 3 {
            simplify::simplify_loop(loop_pts, epsilon)
        } else {
            loop_pts.clone()
        };
        if simplified.len() >= 3 {
            simplified_contours.push(simplified);
        }
    }

    // 3. 全輪郭点を集約 + 制約エッジ作成
    let mut all_points: Vec<(f32, f32)> = Vec::new();
    let mut constraint_edges: Vec<(usize, usize)> = Vec::new();

    for contour in &simplified_contours {
        let base = all_points.len();
        for pt in contour {
            all_points.push(*pt);
        }
        // ループの制約エッジ
        let cn = contour.len();
        for i in 0..cn {
            constraint_edges.push((base + i, base + (i + 1) % cn));
        }
    }

    let _contour_point_count = all_points.len();

    // 4. 内部点サンプリング (Poisson Disk Sampling)
    // パディング後のサイズで処理する
    let seed = (w as u64 * 31337 + h as u64 * 7919 + density as u64 * 104729) | 1;
    let interior_pts = poisson::poisson_disk_sample(
        padded_w,
        padded_h,
        min_spacing,
        &dilated_alpha,
        1, // 膨張済みなのでしきい値1で十分
        &all_points,
        seed,
    );
    all_points.extend_from_slice(&interior_pts);

    // 最低3点必要
    if all_points.len() < 3 {
        let hw = w as f32 * 0.5;
        let hh = h as f32 * 0.5;
        all_points.clear();
        all_points.push((-hw, -hh));
        all_points.push((hw, -hh));
        all_points.push((hw, hh));
        all_points.push((-hw, hh));
        constraint_edges.clear();
    }

    // 5. 制約付き Delaunay 三角形分割
    let tris = delaunay::triangulate(
        &all_points,
        &constraint_edges,
        &dilated_alpha,
        padded_w,
        padded_h,
        1, // 膨張済みなのでしきい値1で十分
    );

    // 結果をC互換メモリに変換
    let vertex_count = all_points.len() as u32;
    let index_count = (tris.len() * 3) as u32;

    let vertices = if vertex_count > 0 {
        let mut verts = Vec::with_capacity(all_points.len() * 2);
        for &(x, y) in &all_points {
            verts.push(x);
            verts.push(y);
        }
        let ptr = verts.as_mut_ptr();
        std::mem::forget(verts);
        ptr
    } else {
        std::ptr::null_mut()
    };

    let indices = if index_count > 0 {
        let mut idxs = Vec::with_capacity(tris.len() * 3);
        for tri in &tris {
            idxs.push(tri[0] as u32);
            idxs.push(tri[1] as u32);
            idxs.push(tri[2] as u32);
        }
        let ptr = idxs.as_mut_ptr();
        std::mem::forget(idxs);
        ptr
    } else {
        std::ptr::null_mut()
    };

    Box::into_raw(Box::new(MeshResult {
        vertices,
        vertex_count,
        indices,
        index_count,
    }))
}

/// MeshResult を解放する。
///
/// Lua側で結果を使い終わった後に呼び出す。
#[no_mangle]
pub extern "C" fn puppet_mesh_free(result: *mut MeshResult) {
    if result.is_null() {
        return;
    }
    unsafe {
        let r = Box::from_raw(result);
        if !r.vertices.is_null() && r.vertex_count > 0 {
            let _ = Vec::from_raw_parts(r.vertices, (r.vertex_count * 2) as usize, (r.vertex_count * 2) as usize);
        }
        if !r.indices.is_null() && r.index_count > 0 {
            let _ = Vec::from_raw_parts(r.indices, r.index_count as usize, r.index_count as usize);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_full_pipeline() {
        // 32x32の円形α画像
        let w: usize = 32;
        let h: usize = 32;
        let hw = w as f32 * 0.5;
        let hh = h as f32 * 0.5;
        let r = 12.0f32;

        let mut rgba = vec![0u8; w * h * 4];
        for y in 0..h {
            for x in 0..w {
                let dx = x as f32 - hw;
                let dy = y as f32 - hh;
                if dx * dx + dy * dy <= r * r {
                    let idx = (y * w + x) * 4;
                    rgba[idx] = 255;     // R
                    rgba[idx + 1] = 255; // G
                    rgba[idx + 2] = 255; // B
                    rgba[idx + 3] = 255; // A
                }
            }
        }

        let result = puppet_mesh_generate(rgba.as_ptr(), w as i32, h as i32, 128, 10, 0);
        assert!(!result.is_null());

        unsafe {
            let r = &*result;
            assert!(r.vertex_count > 0, "頂点が生成されるべき: {}", r.vertex_count);
            assert!(r.index_count > 0, "インデックスが生成されるべき: {}", r.index_count);
            assert_eq!(r.index_count % 3, 0, "インデックス数は3の倍数であるべき");

            // 全インデックスが有効範囲内であることを確認
            for i in 0..r.index_count as usize {
                let idx = *r.indices.add(i);
                assert!(
                    idx < r.vertex_count,
                    "インデックス{}が範囲外: {} >= {}",
                    i,
                    idx,
                    r.vertex_count
                );
            }
        }

        puppet_mesh_free(result);
    }

    #[test]
    fn test_null_input() {
        let result = puppet_mesh_generate(std::ptr::null(), 0, 0, 128, 10, 0);
        assert!(!result.is_null());
        unsafe {
            let r = &*result;
            assert_eq!(r.vertex_count, 0);
            assert_eq!(r.index_count, 0);
        }
        puppet_mesh_free(result);
    }
}
