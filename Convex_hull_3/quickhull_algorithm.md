# CGAL QuickHull 3D 算法深度分析

> 基于 `Convex_hull_3/include/CGAL/convex_hull_3.h` 源码，完整解析 3D 凸包的 QuickHull 实现。

---

## 目录

1. [算法总览](#1-算法总览)
2. [退化情况处理](#2-退化情况处理)
3. [初始四面体的构建](#3-初始四面体的构建)
4. [核心循环：ch_quickhull_3_scan](#4-核心循环ch_quickhull_3_scan)
5. [关键子程序详解](#5-关键子程序详解)
6. [数据结构选择的逻辑](#6-数据结构选择的逻辑)
7. [谓词精度：三层过滤机制](#7-谓词精度三层过滤机制)
8. [复杂度分析](#8-复杂度分析)
9. [与 2D QuickHull（Bykat）的对比](#9-与-2d-quickhullbykat的对比)

---

## 1. 算法总览

QuickHull 3D 的核心思想是**增量扩张**：从一个初始四面体出发，不断找到在当前凸包外部最远的点，将其纳入凸包，同时删除被遮挡的面，添加新的面。

整体流程：

```
输入点集
    │
    ├─ 退化检测（0/1/2点、共线、共面）
    │       └─ 共面 → 调用 convex_hull_2()
    │
    └─ 一般情况
            │
            ├─ 1. 找三个不共线的点 p1, p2, p3
            ├─ 2. 找距离平面(p1,p2,p3)最远的点 p4
            ├─ 3. 构建初始四面体 (p1,p2,p3,p4)
            ├─ 4. 将剩余点分配到四面体各面的"外部集合"
            └─ 5. 主循环：处理所有有外部点的面
                    ├─ 取面 f 的外部集合中最远点 p
                    ├─ 找从 p 可见的所有面（visible set）
                    ├─ 确定可见面的边界轮廓（horizon）
                    ├─ 删除可见面，以 p 为顶点构建新面
                    └─ 将原外部点重新分配到新面
```

---

## 2. 退化情况处理

`convex_hull_3()` 在进入 QuickHull 之前，完整处理所有退化情况（来自 `ch_quickhull_face_graph` 的调用链上层）：

```
点数 = 0          → 返回空
点数 = 1          → 返回 Point_3
点数 = 2          → 返回 Segment_3
点数 = 3 且不共线 → 返回 Triangle_3
所有点共线        → 返回 Segment_3（最远两点）
所有点共面        → 调用 coplanar_3_hull()
                      └─ 选最合适的投影平面（xy/yz/xz）
                      └─ 调用 convex_hull_2()
                      └─ 结果映射回 3D
一般情况          → ch_quickhull_face_graph()
```

共面检测的实现：先找三个不共线的点定义一个平面，然后找距离该平面最远的点。如果最远点也在平面上，则所有点共面。

```cpp
// 来自源码
if (coplanar(*point1_it, *point2_it, *point3_it, *max_it)) {
    coplanar_3_hull(points.begin(), points.end(), ...);
} else {
    // 进入 QuickHull
}
```

---

## 3. 初始四面体的构建

QuickHull 需要一个初始的凸多面体作为起点。CGAL 选择**四面体**（4个顶点，4个三角面）。

**构建过程**（来自 `ch_quickhull_face_graph`）：

```
1. 已有三个不共线的点 p1, p2, p3（在退化检测时找到）
2. 在剩余点中找距离平面(p1,p2,p3)最远的点 p4
   - 用 min_max_element 同时找最大和最小有符号距离
   - 取绝对值更大的那个（可能在平面两侧）
3. 构建四面体的 4 个三角面，保证法向量朝外
```

四面体的面邻接关系（来自源码的硬编码）：

```cpp
Face_handle f0 = tds.create_face(v0, v1, v2);  // 底面
Face_handle f1 = tds.create_face(v3, v1, v0);  // 侧面1
Face_handle f2 = tds.create_face(v3, v2, v1);  // 侧面2
Face_handle f3 = tds.create_face(v3, v0, v2);  // 侧面3

// 邻接关系：f0 的三条边分别与 f2, f3, f1 相邻
f0->set_neighbors(f2, f3, f1);
f1->set_neighbors(f0, f3, f2);
f2->set_neighbors(f0, f1, f3);
f3->set_neighbors(f0, f2, f1);
```

顶点顺序保证每个面的法向量指向四面体外部（逆时针 = 正面朝外）。

**关键设计**：四面体用 `Triangulation_data_structure_2`（TDS_2）存储，而不是 `Surface_mesh` 或 `Polyhedron_3`。原因见第 6 节。

---

## 4. 核心循环：ch_quickhull_3_scan

这是整个算法的心脏，位于 `ch_quickhull_3_scan` 函数中。

### 4.1 数据结构

每个面（`TDS_2::Face`）携带两个额外字段（通过 `Convex_hull_face_base_2` 注入）：

```cpp
std::list<Point_3> points;              // 该面的"外部点集合"
std::list<Face_handle>::iterator it;    // 在 pending_facets 中的位置
```

`pending_facets`：所有外部点集合非空的面的列表，即"还需要处理的面"。

### 4.2 主循环逻辑

```
while pending_facets 非空:
    1. 取 pending_facets 的第一个面 f
    2. 找 f 的外部点集合中距离 f 最远的点 p（farthest_outside_point）
    3. 从 f 出发，BFS 找所有从 p 可见的面（find_visible_set）
       → 输出：visible_set（可见面集合）
       → 输出：border（地平线边，即可见/不可见面的边界）
    4. 收集所有可见面的外部点到 vis_outside_set
    5. 从 pending_facets 中移除所有可见面
    6. 将 border 边排成有序环（通过 map 追踪顶点连接）
    7. 删除可见面，以 p 为顶点，沿 border 环创建新面（star_hole）
    8. 将 vis_outside_set 中的点重新分配到新面（partition_outside_sets）
    9. 将有外部点的新面加入 pending_facets
```

### 4.3 border 边的排序

可见面与不可见面之间的边界（horizon）必须排成有序环，才能用 `star_hole` 填充。CGAL 用一个 `std::map<Vertex_handle, Edge>` 来追踪：

```cpp
// border 中存的是：顶点 v → 以 v 为起点的 horizon 边
// 通过不断查找"当前边终点"作为下一条边的起点，串成有序环
std::map<Vertex_handle, Edge> border;

Edge e = border.begin()->second;
edges.push_back(e);
border.erase(border.begin());
while (!border.empty()) {
    // 找以当前边终点为起点的下一条边
    it = border.find(e.first->vertex(TDS_2::ccw(e.second)));
    e = it->second;
    edges.push_back(e);
    border.erase(it);
}
```

---

## 5. 关键子程序详解

### 5.1 find_visible_set：BFS 找可见面

从已知可见面 `start` 出发，BFS 扩展到所有从点 `p` 可见的面。

**可见性判断**：面 f 对点 p 可见，当且仅当 p 在 f 的法向量正侧（即 `is_on_positive_side` 为 true）。

注意：CGAL 的面法向量约定是**逆时针为正面**，所以"可见"等价于"p 在面的外侧"。

```cpp
// 对每个邻居面 f：
Is_on_positive_side_of_plane_3<Traits> is_on_positive_side(
    traits, f->vertex(0)->point(),
            f->vertex(2)->point(),   // 注意：顶点顺序是 0,2,1（反转法向量）
            f->vertex(1)->point());

if (!is_on_positive_side(point)) {
    // 可见：加入 visible_set，继续 BFS
    visible.push_back(f);
} else {
    // 不可见：这条边是 horizon 边
    border.insert(...);
}
```

**顶点清理**：BFS 过程中同时标记顶点。可见面的顶点如果不在 horizon 上，就从 TDS 中删除（因为它们不再是凸包顶点）。

### 5.2 farthest_outside_point：找最远点

```cpp
// 用 std::max_element + lambda 找距离面最远的点
Outside_set_iterator farthest_it =
    std::max_element(outside_set.begin(), outside_set.end(),
        [&less_dist_to_plane, &plane](const Point_3& p1, const Point_3& p2) {
            return less_dist_to_plane(plane, p1, p2);
        });
```

`Less_signed_distance_to_plane_3` 是 Traits 提供的谓词，比较两点到平面的有符号距离。

### 5.3 partition_outside_sets：重新分配外部点

新面创建后，原来可见面的外部点需要重新分配：

```cpp
for (每个新面 f) {
    Is_on_positive_side_of_plane_3 is_on_positive_side(traits, f的三个顶点);
    for (vis_outside_set 中的每个点 p) {
        if (is_on_positive_side(p)) {
            // p 在 f 的外侧，归属于 f
            f->points.splice(f->points.end(), vis_outside_set, p的迭代器);
        }
    }
}
```

用 `std::list::splice` 移动元素，避免拷贝。一个点只会被分配到第一个包含它的面（因为点一旦被 splice 走就不在 vis_outside_set 里了）。

### 5.4 star_hole：填充 horizon

`tds.star_hole(edges.begin(), edges.end(), visible_set.begin(), visible_set.end())` 是 TDS_2 的内置操作：

- 沿 horizon 边环，以新顶点 p 为中心，创建一圈新三角面
- 复用 visible_set 中的面对象（避免内存分配），多余的删除，不足的新建

```cpp
// 面的数量调整
std::ptrdiff_t diff = visible_set.size() - edges.size();
if (diff < 0) {
    for (int i = 0; i < -diff; i++)
        visible_set.push_back(tds.create_face());  // 新建不足的面
} else {
    for (int i = 0; i < diff; i++) {
        tds.delete_face(visible_set.back());        // 删除多余的面
        visible_set.pop_back();
    }
}
Vertex_handle vh = tds.star_hole(edges.begin(), edges.end(),
                                  visible_set.begin(), visible_set.end());
vh->point() = farthest_pt;
```

---

## 6. 数据结构选择的逻辑

### 为什么用 TDS_2 而不是 Surface_mesh？

`Triangulation_data_structure_2` 是专为三角剖分设计的结构，在 QuickHull 中有几个关键优势：

**1. 面可以携带自定义数据**

通过 `Convex_hull_face_base_2` 模板参数，每个面可以存储 `std::list<Point_3> points`（外部点集合）和 `iterator it`（在 pending_facets 中的位置）。`Surface_mesh` 的 property map 也能做到，但需要额外的 map 查找开销。

**2. `star_hole` 操作**

TDS_2 内置了 `star_hole`，可以高效地以一个新顶点填充一个洞（horizon 围成的区域）。这个操作在 `Surface_mesh` 上需要手动实现。

**3. 面的 `info()` 字段**

TDS_2 的面有 `info()` 字段，用于标记 VISITED/BORDER 状态，BFS 时不需要额外的 `std::set` 或 `std::unordered_set`。

**4. 最终转换**

算法结束后，TDS_2 通过 `copy_face_graph(tds, P)` 转换到用户指定的输出类型（`Surface_mesh`、`Polyhedron_3` 等）。这个转换是 O(n) 的，不影响整体复杂度。

---

## 7. 谓词精度：三层过滤机制

`Is_on_positive_side_of_plane_3` 是整个算法中调用最频繁的谓词，CGAL 对它做了三层精度过滤：

```
第一层：静态过滤（Static Filter）
    ├─ 预计算平面法向量的分量 m10, m20, m21
    ├─ 计算误差界 eps = 5.11e-15 * maxx * maxy * maxz
    ├─ 如果 |det| > eps → 直接返回结果（大多数情况）
    └─ 否则 → 进入第二层

第二层：区间算术（Interval Arithmetic）
    ├─ 用 Interval_nt_advanced 重新计算
    ├─ 如果区间不包含 0 → 确定结果
    └─ 否则 → 进入第三层

第三层：精确算术（Exact Arithmetic）
    └─ 用 Simple_cartesian<Exact_field_selector<double>::Type>
       做精确计算，保证正确性
```

这个三层机制保证了：
- **速度**：99%+ 的情况在第一层就能判断
- **正确性**：退化情况（点恰好在平面上）用精确算术处理

```cpp
// 来自源码的三层结构
int static_res = static_filtered(psx, psy, psz);
if (static_res != STATIC_FILTER_FAILURE)
    return static_res == 1;  // 第一层成功

try {
    Uncertain<Sign> res = sign(scalar_product(...));  // 第二层
    if (is_certain(res)) return (get_certain(res) == POSITIVE);
} catch (Uncertain_conversion_exception&) {}

// 第三层：精确算术
return sign(scalar_product(to_EK(s) - ek_plane_ptr->point,
                           ek_plane_ptr->vector)) == POSITIVE;
```

---

## 8. 复杂度分析

| 情况 | 时间复杂度 | 说明 |
|------|-----------|------|
| 最坏情况 | O(n²) | 每次只添加一个凸包顶点，且每次都要扫描所有剩余点 |
| 平均情况（随机输入） | O(n log n) | 期望每个点只被少数面的外部集合包含 |
| 输出敏感 | O(n log h) | h 为凸包顶点数，理论上可达到，实践中接近 |

**空间复杂度**：O(n)，每个点最多属于一个面的外部集合。

**与 2D 的对比**：2D QuickHull（Bykat）的最坏情况也是 O(n²)，但 3D 的常数更大，因为每次迭代需要 BFS 找可见面集合。

---

## 9. 与 2D QuickHull（Bykat）的对比

| 维度 | 2D Bykat | 3D QuickHull |
|------|---------|-------------|
| 初始结构 | 最左/最右两点（线段） | 四面体（4点4面） |
| 分治方式 | 递归，用栈模拟 | 迭代，用 pending_facets 队列 |
| "外部点"存储 | 用 `std::partition` 原地分割 | 每个面携带 `std::list<Point_3>` |
| 可见性判断 | 点在直线哪侧（`left_turn`） | 点在平面哪侧（`is_on_positive_side`） |
| 边界处理 | 两个端点 | horizon 边环（需要排序） |
| 数据结构 | `std::vector` + 原地操作 | `TDS_2`（三角剖分数据结构） |
| 最终输出 | 有序点序列 | 多面体网格（通过 `copy_face_graph` 转换） |

**核心差异**：3D 的"删除可见面、填充新面"操作在拓扑上比 2D 复杂得多。2D 只需要替换一段弧，3D 需要处理任意形状的 horizon 环，这是引入 TDS_2 和 `star_hole` 的根本原因。

---

## 附录：调用链总览

```
convex_hull_3(first, beyond, P)          ← 用户调用的入口
    │
    └─ ch_quickhull_face_graph(points, p1_it, p2_it, p3_it, P, traits)
            │
            ├─ coplanar_3_hull(...)       ← 共面退化
            │       └─ convex_hull_2()
            │
            └─ 一般情况：
                    ├─ 构建初始四面体（TDS_2）
                    ├─ non_coplanar_quickhull_3(points, tds, traits)
                    │       ├─ 初始点分配到各面
                    │       └─ ch_quickhull_3_scan(tds, pending_facets, traits)
                    │               ├─ farthest_outside_point()
                    │               ├─ find_visible_set()
                    │               ├─ tds.star_hole()
                    │               └─ partition_outside_sets()
                    │
                    └─ copy_face_graph(tds, P)  ← 转换到输出类型
```
