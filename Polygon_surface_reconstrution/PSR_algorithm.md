# CGAL Polygonal Surface Reconstruction Package 算法与实现分析

## 概述

本文档分析 Polygonal Surface Reconstruction 包（又称 PolyFit）的实现细节。该包实现了 Nan & Wonka 2017 论文中的"假设-选择"（hypothesis-and-selection）方法，从点云重建分片平面（piecewise planar）物体的封闭多边形网格。所有实现集中于：

- `include/CGAL/Polygonal_surface_reconstruction.h` — 顶层类及 MIP 问题构建
- `include/CGAL/Polygonal_surface_reconstruction/internal/hypothesis.h` — 候选面生成
- `include/CGAL/Polygonal_surface_reconstruction/internal/compute_confidences.h` — 置信度计算
- `include/CGAL/Polygonal_surface_reconstruction/internal/alpha_shape_mesh.h` — Alpha Shape 覆盖面积估计
- `include/CGAL/Polygonal_surface_reconstruction/internal/point_set_with_planes.h` — 带平面标注的点集

---

## 一、整体流程

算法分为三个阶段，在构造函数和 `reconstruct()` 中分别执行：

```
构造函数：
  1. 平面精化（refine_planes）
  2. 候选面生成（hypothesis_.generate）
  3. 置信度计算（conf.compute）

reconstruct()：
  4. 邻接关系提取（extract_adjacency）
  5. MIP 问题构建与求解
  6. 结果网格重建
```

构造函数完成后，候选面和置信度被缓存在 `candidate_faces_` 中，允许用户以不同权重多次调用 `reconstruct()`，而无需重复执行代价高昂的候选面生成步骤。

---

## 二、候选面生成：`Hypothesis` 类

### 2.1 平面精化：`refine_planes`

在生成候选面之前，先合并近似共面的平面段。判断标准是：如果平面段 s1 上的点有足够多落在平面段 s2 的平面上（在距离阈值 `dist_threshold` 内），则将两者合并。

合并后的平面段使用所有点重新拟合平面，保证平面方程的精度。这一步减少了后续候选面的数量，避免因近似重复平面导致的组合爆炸。

### 2.2 包围盒网格：`construct_bbox_mesh`

为了让平面段在有限范围内相交，算法先构造点集的轴对齐包围盒（AABB），并将其表示为一个六面体网格（`bbox_mesh`）。包围盒的六个面作为"边界平面"参与后续的平面相交计算，确保所有候选面都有界。

### 2.3 代理网格：`construct_proxy_mesh`

对每个平面段，计算其支撑点在平面上的投影，然后用 `CGAL::convex_hull_2` 计算投影点的凸包，得到该平面段在其平面上的初始多边形面。这些初始面被裁剪到包围盒内，形成"代理网格"（proxy mesh）。

### 2.4 平面两两相交：`pairwise_intersection`

这是候选面生成的核心步骤，也是计算量最大的部分。

**顶点计算策略**：每个顶点由三个平面的交点唯一确定。为了数值稳定性，实现预先计算所有平面三元组的交点（`compute_triplet_intersections`），并用三个平面的指针（按地址排序）作为 key 缓存结果（`query_intersection`）。这样同一个几何顶点无论从哪条路径到达，都使用同一个 `Point*`，避免了浮点误差导致的顶点重复。

**切割流程**：对每个平面 P，遍历所有其他平面 Q，用 Q 切割 P 上的所有面：

```
对每个平面 P 上的面 f：
  1. do_intersect(f, Q) 检查 f 是否与 Q 相交
  2. compute_intersections 计算交点（已有顶点或边上的新顶点）
  3. split_edge 在边上插入新顶点（使用 Euler::split_edge）
  4. cut 将面沿交线分割为两个子面
```

`cut` 函数内部使用 `Euler::split_face` 操作，在半边数据结构上直接修改拓扑，保证网格的一致性。

**结果**：每个平面段被分割为若干多边形候选面，每个候选面携带其所属平面段的信息（通过 `f:supp_segment` 属性存储）。

### 2.5 邻接关系提取：`extract_adjacency`

邻接关系描述了哪些候选面共享同一条边（即由同一对平面的交线产生）。每个 `Intersection` 结构存储共享该边的所有半边，以及边的两个端点指针。

关键约束：只有恰好有 4 个面共享一条边的情况（`fan.size() == 4`）才会产生 MIP 变量。边界边（`fan.size() < 4`）被强制设为 0（不允许开放边界，保证封闭性）。

---

## 三、置信度计算：`Candidate_confidences` 类

置信度计算为每个候选面计算三个属性，存储为网格的 Property Map：

### 3.1 支撑点数：`f:num_supporting_points`

`supporting_points` 函数找出属于该面所在平面段的所有点，然后过滤出投影落在面内部的点：

```cpp
// 将点投影到平面，检查是否在面的 2D 多边形内
Point2 proj = plane.to_2d(point_3d);
if (polygon_2d.bounded_side(proj) == ON_BOUNDED_SIDE)
    indices.push_back(point_index);
```

这里使用 Polygon 包的 `bounded_side_2` 做点包含测试，体现了包之间的复用。

### 3.2 面积：`f:face_area`

对多边形面做扇形三角剖分（以第一个顶点为中心），累加各三角形面积：

```cpp
result += triangle_area(p, q, r);  // 叉积模长的一半
```

### 3.3 覆盖面积：`f:covered_area`

这是置信度计算中最复杂的部分，用于估计候选面中被点云实际覆盖的面积（即 `area(M_i^α)`）。

**方法**：Alpha Shape。将支撑点投影到候选面所在平面，构造 2D Alpha Shape，提取 Alpha Shape 的内部三角形，将其映射回 3D 空间，计算总面积。

Alpha 值的选取：使用点云的平均间距（`CGAL::compute_average_spacing`）的两倍作为 Alpha 值。这是一个经验参数，保证 Alpha Shape 能覆盖点云密集区域而不过度扩张到稀疏区域。

`Alpha_shape_mesh` 类封装了这一过程：
1. 将 3D 点投影到平面（`plane.to_2d`）
2. 构造带层次结构的 Delaunay 三角化（`Triangulation_hierarchy_2`，加速点定位）
3. 构造 Alpha Shape（`Alpha_shape_2`）
4. 提取 `INTERIOR` 类型的三角形，映射回 3D 坐标

---

## 四、MIP 问题构建：`reconstruct()` 方法

### 4.1 变量设计

MIP 问题有三类二值变量（`Variable::BINARY`）：

| 变量范围 | 含义 |
|---|---|
| `x[0..N-1]` | 第 i 个候选面是否被选中（`x_i ∈ {0,1}`） |
| `x[N..N+E-1]` | 第 j 条内部边是否被使用（`edge_usage`） |
| `x[N+E..N+2E-1]` | 第 j 条内部边是否为锐利边（`edge_sharp`） |

只有 `fan.size() == 4` 的边才产生后两类变量（边界边不产生，因为它们被强制为 0）。

### 4.2 目标函数

目标函数是三个能量项的加权和，实现时对系数做了归一化处理以平衡量纲：

```cpp
// 数据拟合项系数（直接使用权重）
double coeff_data_fitting = wt_fitting;

// 覆盖项系数（归一化到点数/包围盒面积）
double coeff_coverage = total_points * wt_coverage / box_area;

// 复杂度项系数（归一化到点数/边数）
double coeff_complexity = total_points * wt_complexity / double(adjacency.size());
```

对每个候选面 f，目标函数贡献为：

```
-coeff_data_fitting * num_supporting_points[f]   // 数据拟合（负号：支撑点多的面倾向被选）
+ coeff_coverage * (face_area[f] - covered_area[f])  // 覆盖项（未覆盖面积越小越好）
```

对每条内部边，复杂度项贡献为：

```
+coeff_complexity * x[edge_sharp_idx]  // 锐利边越少越好
```

### 4.3 硬约束

**约束一：流形约束**（每条边恰好被 0 或 2 个面使用）

```
Σ x[f_j] - 2 * x[edge_usage_idx] = 0   （对每条内部边）
Σ x[f_j] = 0                            （对每条边界边，强制不选）
```

这保证了最终网格是 2-流形：每条边恰好连接两个面（或不被使用）。

**约束二：锐利边约束**（若边被标记为锐利，则该边必须被使用）

```
x[edge_usage_idx] - x[edge_sharp_idx] >= 0
```

**约束三：锐利边定义约束**（若两个来自不同平面的面都被选中且共享该边，则该边是锐利的）

使用大 M 法（Big-M method）线性化"若 A 且 B 则 C"的逻辑：

```
x[edge_sharp_idx] - M*x[fid1] - M*x[fid2] - M*x[edge_usage_idx] >= 1 - 3M
```

其中 M = 1.0。当 `x[fid1] = x[fid2] = x[edge_usage_idx] = 1` 时，约束变为 `x[edge_sharp_idx] >= 1`，即强制该边为锐利边。

### 4.4 求解与结果提取

调用 `MixedIntegerProgramTraits::solve()` 求解 MIP 问题。求解成功后：

1. 读取解向量 `X`，将 `round(X[f_idx]) == 0` 的面删除（`Euler::remove_face`）。
2. 标记锐利边（`e:sharp_edges` 属性）。
3. 将内部 `Surface_mesh` 转换为输出网格：
   - `polygon_mesh_to_polygon_soup`：转为点+面索引格式
   - `merge_duplicate_points_in_polygon_soup`：合并重复顶点
   - `orient_polygon_soup`：统一面朝向
   - `polygon_soup_to_polygon_mesh`：转回网格格式

最后两步使用 PMP 包的工具函数，体现了包之间的复用。

---

## 五、数值稳定性处理

### 5.1 平面三元组交点缓存

平面三元组的交点用三个平面指针（按地址升序排列）作为 key 缓存：

```cpp
// 三个平面指针排序后作为 map key
std::tuple<const Plane*, const Plane*, const Plane*> key =
    sorted_triple(plane1, plane2, plane3);
```

这确保同一几何顶点在整个算法中只计算一次，所有引用该顶点的面都使用同一个 `Point*`，避免了浮点误差导致的顶点不一致。

### 5.2 浮点/精确数值的 CDT Tag 自动选择

与 Polygon Repair 包类似，`Hypothesis` 中的三角化操作根据 `Kernel::FT` 是否为浮点类型自动选择相交处理策略，保证在使用 EPICK 时也能正确处理平面相交的退化情况。

---

## 六、性能特征

| 阶段 | 复杂度 | 瓶颈 |
|---|---|---|
| 平面精化 | O(P²) | P 为平面数，通常很小 |
| 候选面生成 | O(P² · F) | F 为每个平面上的面数，受相交次数影响 |
| 置信度计算 | O(N · n) | N 为候选面数，n 为点数 |
| MIP 求解 | 指数级（最坏） | 实际取决于求解器和问题规模 |

MIP 求解是整个流程的性能瓶颈。变量数为 `N + 2E`（N 为候选面数，E 为内部边数），约束数约为 `3E`。对于超过 5000 个候选面的问题，建议使用 Gurobi 或 CPLEX 等商业求解器，GLPK 只适合小规模问题（约 1000 个变量以内）。
