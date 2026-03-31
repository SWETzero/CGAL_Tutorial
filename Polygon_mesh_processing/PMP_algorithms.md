# CGAL Polygon Mesh Processing Package 算法与实现分析

## 概述

本文档分析 PMP 核心包中各算法的实现细节，涵盖法向量计算、度量计算、自相交检测、曲率计算、连通分量分析等。所有实现位于：

- `include/CGAL/Polygon_mesh_processing/compute_normal.h` — 法向量计算
- `include/CGAL/Polygon_mesh_processing/measure.h` — 几何度量（边长、面积、体积）
- `include/CGAL/Polygon_mesh_processing/self_intersections.h` — 自相交检测
- `include/CGAL/Polygon_mesh_processing/curvature.h` — 离散曲率
- `include/CGAL/Polygon_mesh_processing/connected_components.h` — 连通分量
- `include/CGAL/Polygon_mesh_processing/distance.h` — 网格间距离

---

## 一、法向量计算：`compute_normal.h`

### 1.1 面法向量：`compute_face_normal`

**算法**：Newell 方法（适用于任意多边形面，不仅限于三角形）。

对于三角形面，退化为叉积：

```cpp
// triangle_normal 内部实现
Vector_3 n = cross_product(p1→p2, p1→p0);
// 返回值已乘以 1/2，即三角形有向面积向量
return scale(n, FT(1)/FT(2));
```

对于一般多边形面，遍历所有顶点，累加相邻顶点对的叉积贡献（Newell 公式）：

```cpp
// 对面的每条半边 h
Vector_3 n = cross_product(p_prev→p_curr, p_prev→p_next);
normal = normal + n;
```

最终结果通过 `normalize` 函数单位化。`normalize` 使用 `CGAL::approximate_sqrt` 计算模长，并在模长为零时跳过归一化（避免 NaN）。

### 1.2 顶点法向量：`compute_vertex_normal`

**算法**：面积加权平均（Area-Weighted Average）。

遍历顶点的所有相邻面，将每个面的**未归一化**法向量（即面积向量）累加，最后归一化：

```cpp
// 遍历顶点 v 的所有相邻面
for each face f around v:
    normal += triangle_normal(p0, p1, p2);  // 未归一化，模长 = 面积
// 最后归一化
normalize(normal);
```

这等价于按面积加权：面积大的面对顶点法向量贡献更大。这是最常用的顶点法向量估计方法，在光滑曲面上效果良好。

**边界顶点处理**：对于边界顶点，只累加有效面（非 null face）的法向量贡献，边界半边对应的 null face 被跳过。

### 1.3 批量计算：`compute_face_normals` / `compute_vertex_normals`

批量版本遍历网格的所有面/顶点，将结果写入 Property Map。实现上直接调用单元素版本，无特殊优化。

---

## 二、几何度量：`measure.h`

### 2.1 边长：`edge_length`

直接计算两端点的欧氏距离：

```cpp
return approximate_sqrt(squared_distance(source_point, target_point));
```

使用 `approximate_sqrt` 而非精确平方根，因为边长本身就是一个近似量（除非使用精确构造 Kernel）。

### 2.2 面积：`face_area`

**三角形面**：利用叉积模长的一半：

```cpp
// 叉积模长 = 平行四边形面积 = 2 × 三角形面积
return approximate_sqrt(squared_area) / 2;
```

**一般多边形面**：三角剖分后累加各三角形面积（以第一个顶点为扇形中心）：

```cpp
// 与 Polygon 包的 polygon_area_2 算法相同，但在 3D 中
for each triangle (v0, vi, vi+1) in fan triangulation:
    area += face_area(triangle);
```

### 2.3 体积：`volume`

**前提**：网格必须是封闭的（无边界）。

**算法**：散度定理（Divergence Theorem）的离散版本。对每个三角形面，计算其对总体积的有符号贡献：

```cpp
// 对每个三角形 (p0, p1, p2)
// 贡献 = (1/6) * dot(p0, cross(p1, p2))
// 这等价于从原点到三角形的有向四面体体积
volume += dot(p0, cross(p1, p2)) / 6;
```

结果的符号取决于网格的朝向（外法向量朝外为正）。这与 Polygon 包中 Shoelace 公式的三维推广在数学上是等价的。

### 2.4 面匹配：`match_faces`

这是 `measure.h` 中一个不太显眼但很实用的函数，用于判断两个网格的面是否几何上重合（即使拓扑不同）。

**算法**：
1. 对每个面，计算其顶点 id 的规范化排列（最小 id 在前，通过 `rearrange_face_ids` 实现）。
2. 用规范化排列作为 key，建立哈希表。
3. 在另一个网格中查找匹配。

```cpp
// rearrange_face_ids: 将最小元素旋转到首位
auto min_elem = std::min_element(ids.begin(), ids.end());
std::rotate(ids.begin(), min_elem, ids.end());
```

---

## 三、自相交检测：`self_intersections.h`

### 3.1 整体架构

自相交检测是 PMP 中实现最复杂的算法之一，支持三种输入：
- Triangle Mesh（半边网格）
- Triangle Soup（无拓扑三角形集合）
- 两个网格之间的相交检测

通过 `Triangle_mesh_and_triangle_soup_wrapper` 适配器统一处理两种输入格式。

### 3.2 算法：Box Intersection + 精确三角形相交测试

**第一阶段（粗筛）**：用 `CGAL::box_intersection_d` 找出所有 AABB（轴对齐包围盒）相交的面对。

```cpp
// 为每个面构造 AABB
Box box(face_bbox, face_descriptor);
// 用 box_intersection_d 找出所有相交的 box 对
CGAL::box_intersection_d(boxes.begin(), boxes.end(),
                          boxes.begin(), boxes.end(),
                          output_iterator);
```

**第二阶段（精确测试）**：对每对 AABB 相交的面，排除共享顶点/边的相邻面对，然后做精确的三角形-三角形相交测试：

```cpp
// 排除相邻面（共享边或顶点）
if (faces_have_a_shared_edge(f, g, vh, tm)) return;

// 精确相交测试
if (do_intersect(triangle_f, triangle_g))
    output(f, g);
```

### 3.3 相邻面的排除逻辑

`faces_have_a_shared_edge` 函数检查两个面是否共享一条边（两个顶点）。实现上遍历第一个面的三条半边，检查其对面是否等于第二个面：

```cpp
halfedge_descriptor h = halfedge(f, tm);
for (int i = 0; i < 3; ++i) {
    if (face(opposite(h, tm), tm) == g) return true;
    h = next(h, tm);
}
```

仅共享一个顶点的面不被排除，因为它们在几何上可能真的相交。

### 3.4 并行化

当链接 TBB 时，面对的遍历可以并行化：

```cpp
#ifdef CGAL_LINKED_WITH_TBB
tbb::parallel_for(tbb::blocked_range<std::size_t>(0, face_pairs.size()),
    [&](const tbb::blocked_range<std::size_t>& r) {
        for (std::size_t i = r.begin(); i != r.end(); ++i)
            process_pair(face_pairs[i]);
    });
#endif
```

结果收集使用 `tbb::concurrent_vector` 保证线程安全。

---

## 四、曲率计算：`curvature.h` 与 `interpolated_corrected_curvatures.h`

### 4.1 离散曲率（`curvature.h`）

提供基于 cotangent 权重的离散 Laplace-Beltrami 算子计算均值曲率，以及基于角亏量（angle defect）的高斯曲率。

**均值曲率**（Mean Curvature）：

```
H(v) = (1 / 2A) * |Σ (cot α_ij + cot β_ij) * (v - v_j)|
```

其中 `A` 是顶点的混合面积（Voronoi 面积或钝角三角形的修正面积），`α_ij`、`β_ij` 是边 `(v, v_j)` 对面的两个角。

**高斯曲率**（Gaussian Curvature）：

```
K(v) = (2π - Σ θ_i) / A
```

其中 `θ_i` 是顶点 `v` 处各相邻三角形的内角，`A` 是混合面积。

### 4.2 插值修正曲率（`interpolated_corrected_curvatures.h`）

这是 PMP 中较新的曲率估计方法，基于 Lachaud et al. 2020 的论文，通过对法向量场的插值修正来提高曲率估计的精度。

核心思想：在每个三角形内，用顶点法向量的线性插值构造一个修正的法向量场，然后对该场求散度/旋度得到曲率。这种方法对网格质量的鲁棒性更强。

---

## 五、连通分量：`connected_components.h`

### 5.1 面连通分量：`connected_components`

**算法**：BFS/DFS 遍历，通过半边的 `opposite` 操作跨越面边界。

```cpp
// 从种子面出发，BFS 扩展
std::queue<face_descriptor> queue;
queue.push(seed_face);
while (!queue.empty()) {
    face_descriptor f = queue.front(); queue.pop();
    // 遍历 f 的所有相邻面
    for each halfedge h of f:
        face_descriptor neighbor = face(opposite(h, tm), tm);
        if (neighbor != null_face && !visited[neighbor]):
            component_id[neighbor] = current_id;
            queue.push(neighbor);
}
```

**约束边**：可以通过 `edge_is_constrained_map` 指定不可跨越的约束边，使得约束边两侧的面属于不同连通分量。这在处理特征线（sharp edges）时非常有用。

### 5.2 体积连通分量：`volume_connected_components`

对于封闭网格，还可以计算体积连通分量（即被封闭曲面包围的体积区域）。这需要先确定每个面的朝向，然后追踪哪些体积区域是连通的。

---

## 六、距离计算：`distance.h`

### 6.1 点到网格距离：`max_distance_to_triangle_mesh`

使用 AABB 树加速最近点查询：

```cpp
AABB_tree tree(faces(tm).first, faces(tm).second, tm);
tree.accelerate_distance_queries();

// 对每个查询点
for each point p in point_range:
    max_dist = max(max_dist, tree.squared_distance(p));
```

### 6.2 Hausdorff 距离：`bounded_error_Hausdorff_distance`

PMP 实现了一种有界误差的 Hausdorff 距离近似算法（Attene 2020），保证结果在真实 Hausdorff 距离的 `error_bound` 范围内：

```
|approx_H - true_H| ≤ error_bound
```

**算法核心**：
1. 对每个面，计算其到另一网格的最大距离的下界和上界。
2. 如果上界 - 下界 ≤ error_bound，停止细化。
3. 否则，细分该面并递归处理。

这是一种自适应细化策略，在需要高精度的区域自动加密采样，在平坦区域减少计算量。

---

## 七、内部实现模式

### 7.1 `internal` 命名空间

与 Polygon 包类似，PMP 将所有内部实现放在 `CGAL::Polygon_mesh_processing::internal` 命名空间下，不作为公开 API。用户代码不应直接调用这些函数。

### 7.2 `choose_parameter` 模式

Named Parameters 的提取统一使用 `choose_parameter` 函数：

```cpp
// 从 named parameters 中提取 vertex_point_map，若未提供则使用默认值
auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, pmesh));
```

这个模式在 PMP 的每个函数入口处几乎都会出现，是理解 PMP 源码的关键模式。

### 7.3 `GetGeomTraits` 萃取器

PMP 通过 `GetGeomTraits` 从 Named Parameters 或网格类型中自动推导几何 Traits：

```cpp
typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
```

用户通常不需要显式指定 Traits，PMP 会从网格的 `vertex_point_map` 的值类型自动推导出对应的 Kernel。
