# CGAL Polygon Repair Package 算法与实现分析

## 概述

本文档分析 Polygon Repair 包中各核心算法的实现细节，涵盖奇偶规则修复、非零绕数修复、布尔运算、有效性验证等。所有实现位于：

- `include/CGAL/Polygon_repair/repair.h` — 顶层 `repair`/`is_valid`/`join`/`intersect` 函数及 `Polygon_repair` 类
- `include/CGAL/Polygon_repair/Winding.h` — Non-zero 规则的绕数计算后端
- `include/CGAL/Polygon_repair/Boolean.h` — Union/Intersection 规则的布尔运算后端
- `include/CGAL/Polygon_repair/internal/Triangulation_with_even_odd_constraints_2.h` — Even-odd 规则的特殊约束插入

---

## 一、Even-Odd 规则修复

### 1.1 约束插入：`even_odd_insert_constraint`

Even-odd 规则的核心在于约束的插入方式，实现于 `Triangulation_with_even_odd_constraints_2`。

**关键思想**：如果一条边已经存在于三角化中，则**删除**它而不是重复插入。这实现了奇偶抵消：同一条边出现偶数次等于不存在，出现奇数次等于存在一次。

```cpp
void even_odd_insert_constraint(Vertex_handle va, Vertex_handle vb) {
    // 如果 [va, vc] 是已有边（vc 是 vb 方向上的第一个顶点）
    if (includes_edge(va, vb, vc, incident_face, opposite_vertex)) {
        if (is_constrained(Edge(incident_face, opposite_vertex)))
            remove_constrained_edge(incident_face, opposite_vertex); // 已有 → 删除
        else
            mark_constraint(incident_face, opposite_vertex);          // 没有 → 添加
        if (vc != vb) even_odd_insert_constraint(vc, vb); // 递归处理剩余部分
        return;
    }
    // 如果 [va, vb] 与已有约束边相交，先找到交点再递归
    if (find_intersected_faces(..., intersection)) {
        even_odd_insert_constraint(va, intersection);
        even_odd_insert_constraint(intersection, vb);
        return;
    }
    // 否则正常插入
    triangulate_hole(intersected_faces, ...);
}
```

这个递归实现处理了三种情况：边已存在、边与已有约束相交、边可以直接插入。

### 1.2 唯一边集合的预处理

在调用 `even_odd_insert_constraint` 之前，`add_to_triangulation_even_odd` 先对输入边做预处理，用 `unique_edges`（一个 `unordered_set` 或 `set`）实现奇偶抵消：

```cpp
// 对每条边，规范化方向（小端点在前）后尝试插入集合
// 如果已存在则删除（偶数次抵消），不存在则插入（奇数次保留）
std::pair<Point_2, Point_2> pair = (edge.source() < edge.target()) ?
    make_pair(source, target) : make_pair(target, source);
auto inserted = unique_edges.insert(pair);
if (!inserted.second) unique_edges.erase(inserted.first); // 已存在 → 删除
```

这是一个在插入三角化之前的快速预过滤，避免将重复边插入 CDT。

### 1.3 标签传播：`label_triangulation_even_odd`

标签传播分两步：

**步骤一：共线顶点简化**。遍历所有顶点，如果一个顶点恰好有两条约束边且三点共线，则删除该顶点并用一条长约束边替代。这消除了因边分割产生的中间顶点，避免标签传播时的顺序依赖问题。

**步骤二：BFS 标签传播**。从无穷面（标签 -1，代表外部）出发，遇到约束边时不穿越，而是将对面的三角形加入待检查队列。队列中的三角形根据"是从外部区域还是内部区域穿越约束边到达的"来决定标签：

```cpp
// 从外部（label < 0）穿越约束边 → 新区域是内部（label = +N）
if (to_check_added_by.front() < 0)
    label_region(t, face, ++number_of_polygons, ...);
// 从内部（label > 0）穿越约束边 → 新区域是外部/洞（label = -M）
else
    label_region(t, face, -(++number_of_holes + 1), ...);
```

`label_region` 本身是一个 BFS，在不穿越约束边的情况下将连通区域内所有三角形标记为同一标签。

### 1.4 多边形重建：`reconstruct_multipolygon`

重建阶段遍历所有有限面，找到标签为正（内部）且未处理的面，对每个面检查其三条边是否与相邻面标签不同（即边界边），从边界边出发调用 `reconstruct_ring` 重建环。

`reconstruct_ring` 的核心是沿边界行走：

```cpp
do {
    ring.push_back(pivot_vertex->point());
    // 绕 pivot_vertex 旋转，找到下一个同标签的面
    Face_circulator fc = t.incident_faces(pivot_vertex, current_face);
    do { ++fc; } while (fc->label() != current_face->label());
    current_face = fc;
    current_opposite_vertex = fc->cw(fc->index(pivot_vertex));
} while (current_face != face_adjacent_to_boundary || ...);
```

重建后通过 `polygon.orientation()` 判断是外边界（逆时针）还是洞（顺时针），分别存入对应的容器。最终用 `Polygon_with_holes_less` 比较器将结果存入 `std::set` 实现字典序排序。

---

## 二、Non-Zero 规则修复：`Winding` 类

### 2.1 绕数计算：`label`

Non-zero 规则使用 `Constrained_triangulation_plus_2`（CDT+），它能追踪每条约束边的方向和重数（通过 `contexts` 接口）。

绕数计算从无穷面（绕数 = 0）出发，BFS 传播：

```cpp
// 穿越约束边时，计算方向增量 delta
for (Context c : cdt.contexts(u, v)) {
    if (*c.current() == u && *next(c.current()) == v)
        ++delta;  // 边方向与穿越方向相同 → +1
    else if (*c.current() == v && *next(c.current()) == u)
        --delta;  // 边方向与穿越方向相反 → -1
}
border.push_back({Edge(n, n->index(fh)), index + delta});
```

`delta` 的计算通过 `contexts` 遍历所有经过该边的约束（可能有多条重叠约束），累加方向贡献。这正确处理了部分重叠边的情况。

### 2.2 域标签：`label_domains`

绕数计算完成后，`label_domains` 将绕数非零的三角形聚合为连通域，每个连通域分配一个唯一标签（用于后续重建时区分不同多边形）：

```cpp
for (auto const face: cdt.all_face_handles()) {
    if (face->info().wind != 0 && face->info().label == 0)
        label_domain(face, label++, ...);
}
```

`label_domain` 是一个 BFS，在绕数非零的三角形之间传播，不穿越绕数为零的三角形（外部）。

### 2.3 与 Even-odd 的实现差异

| 方面 | Even-odd（`Polygon_repair`） | Non-zero（`Winding`） |
|---|---|---|
| CDT 类型 | `Triangulation_with_even_odd_constraints_2` | `Constrained_triangulation_plus_2` |
| 约束插入 | 奇偶抵消插入 | 普通插入（保留方向信息） |
| 面属性 | `label`（正/负整数） | `wind`（绕数）+ `label`（域 ID） |
| 内部判断 | `label > 0` | `wind != 0` |

---

## 三、Union/Intersection 规则：`Boolean` 类

### 3.1 覆盖层数计算：`mark_domains`

`Boolean` 类的核心是计算每个三角形被多少个输入多边形覆盖（`layers` 字段）。算法与 `Winding::label` 类似，但语义不同：

```cpp
// 穿越约束边时，delta 表示覆盖层数的变化
// 进入多边形内部 → +1，离开 → -1
mark_domains(n, fh->info().layers + delta, border);
```

对于两个输入多边形：
- `layers = 0`：不在任何多边形内（外部）
- `layers = 1`：在恰好一个多边形内
- `layers = 2`：在两个多边形的重叠区域内

### 3.2 规则通过谓词函数对象实现

Union 和 Intersection 规则通过传入不同的谓词函数对象来选择"内部"三角形：

```cpp
// Union 规则：layers > 0（在至少一个多边形内）
struct Larger_than_zero {
    bool operator()(int i) const { return i > 0; }
};

// Intersection 规则：layers == N（在所有 N 个多边形内）
struct Equal {
    int val;
    bool operator()(int i) const { return i == val; }
};

// 通过 operator() 触发重建
return bops(Larger_than_zero{});  // Union
return bops(Equal{N});            // Intersection
```

`Boolean::operator()(fct)` 先调用 `label_domains(fct)` 将满足谓词的三角形聚合为连通域，再调用 `reconstruct_ring` 重建边界。这与 `Polygon_repair` 的重建逻辑几乎相同，但代码目前是复制的（源码中有 `@todo taken from Polygon_repair should be factorized` 注释）。

---

## 四、有效性验证：`is_valid`

### 4.1 `Polygon_2` 的验证

最简单的验证，三个条件：
1. 顶点数 ≥ 3
2. 无退化边（`source == target`）
3. 多边形是简单的（调用 Polygon 包的 `is_simple()`）

### 4.2 `Polygon_with_holes_2` 的验证

使用 `Constrained_triangulation_plus_2` 配合 `No_constraint_intersection_tag`（约束相交时抛出异常）进行验证：

**阶段一：外边界验证**。将外边界插入 CDT，若有自相交则捕获异常返回 false。

**阶段二：洞的位置验证**。将所有洞插入同一 CDT，然后从无穷面出发做标签传播，但**只穿越外边界约束**（通过 `Outer_hull_constraint_check` 函数对象过滤）。传播结束后，检查每条洞约束边的两侧是否都在外边界内部（标签 = 1）：

```cpp
// 洞边界的两侧都应该在外边界内部
if (fh->label() != 1 || nfh->label() != 1) {
    std::cerr << "Invalid: outward hole edge" << std::endl;
    return false;
}
```

**阶段三：内部连通性验证**。重置标签，这次穿越所有约束边（包括洞边界），验证内部只有一个连通区域（`regions == 1`）。

### 4.3 `Multipolygon_with_holes_2` 的验证

在验证每个子多边形有效后，将所有多边形的约束插入同一 CDT，逐一验证每个多边形的约束边是否与其他多边形的内部不重叠：

```cpp
// 对第 i 个多边形的约束边，检查其两侧是否有其他多边形的内部
if (hits == 1) { // 是第 i 个多边形的边界
    if (fh->label() > 0 && nfh->label() > 0) {
        // 两侧都是内部 → 该边界完全在另一个多边形内部 → 无效
        return false;
    }
}
```

---

## 五、`Triangulation_face_base_with_repair_info_2`

这是 Even-odd 后端专用的面基类，在标准 `Constrained_triangulation_face_base_2` 基础上增加了两个字段：

```cpp
int label_;      // 区域标签：正数=内部多边形，负数=外部/洞，0=未标记
bool processed_; // 是否已加入处理队列（防止重复入队）
```

`label()` 和 `processed()` 提供读写访问。这两个字段在标签传播和多边形重建两个阶段都被复用（重建前会重置 `processed_` 为 false）。

---

## 六、性能特征

修复算法的时间复杂度主要由两部分决定：

- **CDT 构建**：O(n log n + k log n)，其中 n 为顶点数，k 为约束边相交数。
- **标签传播与重建**：O(n + k)，线性于三角化的面数。

文档中的性能数据（10 万顶点、约 300 个洞的多边形在 0.65 秒内完成）表明，对于"相交数有限"的实际数据集，算法性能良好。最坏情况下（所有边两两相交）k = O(n²)，性能会显著下降。
