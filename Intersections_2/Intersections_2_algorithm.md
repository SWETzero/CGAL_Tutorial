# CGAL Intersections_2 Package 算法与实现分析

## 概述

本文档分析 Intersections_2 包中各核心算法的实现细节，涵盖 2D 几何对象之间的相交检测（`do_intersect`）与相交计算（`intersection`）。所有算法实现位于以下文件中：

- `include/CGAL/intersection_2.h` — 公开 API 入口，宏展开分发
- `include/CGAL/Intersection_traits_2.h` — 返回类型 traits，定义各对象对的交集类型
- `include/CGAL/Intersections_2/XxxType_YyyType.h` — 各对象对的具体实现
- `include/CGAL/Intersections_2/internal/` — 复杂算法的内部实现

---

## 一、架构设计：分发机制与返回类型

### 宏分发

每个对象对的头文件末尾都有两个宏调用：

```cpp
CGAL_INTERSECTION_FUNCTION(Segment_2, Line_2, 2)
CGAL_DO_INTERSECT_FUNCTION(Segment_2, Line_2, 2)
```

这两个宏展开后生成全局函数 `intersection(seg, line, k)` 和 `do_intersect(seg, line, k)`，将调用转发到 `Intersections::internal` 命名空间中的具体实现。对于自身与自身相交（如 `Segment_2` 与 `Segment_2`），使用 `_SELF` 变体。

### `Intersection_traits_2`

`Intersection_traits_2<K, A, B>` 通过模板特化定义了每对类型的交集结果类型。例如：

```cpp
// Segment_2 ∩ Segment_2 → Point_2 或 Segment_2
template<typename K>
struct Intersection_traits<K, typename K::Segment_2, typename K::Segment_2> {
  typedef typename boost::variant<typename K::Point_2, typename K::Segment_2> variant_type;
  typedef boost::optional<variant_type> result_type;
};
```

返回值是 `std::optional<std::variant<...>>`，调用方通过 `std::get<>` 或 `boost::get<>` 提取具体类型。

---

## 二、线段与线段相交：`Segment_2_Segment_2`

**文件**：`Intersections_2/Segment_2_Segment_2.h`

**返回类型**：`Point_2` 或 `Segment_2`（共线重叠时）

### 核心数据结构：`S2S2_inter_info`

算法内部用 `S2S2_inter_info` 结构体传递相交的组合信息，避免重复计算：

```cpp
struct S2S2_inter_info {
  bool inter = false;   // 是否相交
  bool dim = 0;         // 0=点, 1=线段
  std::array<int, 2> pt_ids = {-1,-1}; // 交点是哪个输入端点（-1表示内部交点）
  int config;           // 端点字典序配置，用于确定性地计算交点
};
```

`pt_ids` 中 0/1 表示 seg1 的端点，2/3 表示 seg2 的端点。这个设计使得当交点恰好是某个输入端点时，可以直接返回该端点，避免浮点计算误差。

### 算法：`do_intersect_with_info`

**复杂度**：O(1)，仅使用谓词，无浮点除法。

**核心思想**：先对两条线段的端点按字典序排序（`Less_xy_2`），将所有情况归约为两种拓扑构型，再用 `orientation_2` 谓词判断。

**步骤**：

1. 将 seg1 的端点按字典序排为 A1 < A2，seg2 的端点排为 B1 < B2。
2. 用 bbox 快速过滤：若 A2 < B1 或 B2 < A1，则不相交。
3. 按 A1 与 B1 的字典序关系分三大类（SMALLER / EQUAL / LARGER），每类再按端点的相对位置分子情况：

**情况一：A1 < B1（crossing 构型）**，即 A1 < B1 < A2 < B2：

```
A1 -------- A2
       B1 -------- B2
```

调用 `seg_seg_do_intersect_crossing`，检查 B1 在 A1A2 的哪侧，以及 A2 在 B1B2 的哪侧：

```cpp
// 若 B1 在 A1A2 左侧（LEFT_TURN），则 A2 必须在 B1B2 右侧才相交
switch (orientation(A1, A2, B1)) {
  case LEFT_TURN:
    return orientation(B1, B2, A2) == RIGHT_TURN;
  case RIGHT_TURN:
    return orientation(B1, B2, A2) == LEFT_TURN;
  case COLLINEAR:
    // B1 在 A1A2 上，交点就是 B1
    return S2S2_inter_info(B1_id);
}
```

**情况二：A1 < B1 < B2 < A2（contained 构型）**，seg2 完全被 seg1 包含：

```
A1 -------------- A2
      B1 ---- B2
```

调用 `seg_seg_do_intersect_contained`，检查 B1 和 B2 分别在 A1A2 的哪侧。

### 交点计算：`s2s2_alpha`

当两线段在内部相交时，用参数化方法计算交点。参数 α 表示交点在 seg1 上的位置（α=0 为 A1，α=1 为 A2）：

```
α = [(B1 - A1) × (B2 - B1)] / [(A2 - A1) × (B2 - B1)]
```

其中 × 表示 2D 叉积（行列式）。对 `double` 类型使用 `std::fma` 提高精度：

```cpp
double val = std::fma(lx, s2_dy, -ly*s2_dx) / std::fma(s1_dx, s2_dy, -s1_dy*s2_dx);
```

交点通过 `construct_barycentric_2_object()(A1, alpha, A2)` 构造，即 `(1-α)*A1 + α*A2`。

**特殊优化**：对于轴对齐线段（水平或垂直），直接用坐标构造交点，避免参数化计算：

```cpp
// 一条垂直（x相同），一条水平（y相同）
if (pts[0].x() == pts[1].x() && pts[2].y() == pts[3].y())
  _intersection_point = Point_2(pts[0].x(), pts[2].y());
```

### `config` 字段的作用

`config` 是 0-7 的整数，编码了四个端点的字典序关系。它决定了在计算内部交点时，应该用哪条线段作为"主线段"（参数化的基准），以及端点的顺序。这保证了对于相同的两条线段，无论以何种顺序传入，计算出的交点坐标是确定性的（deterministic）。

---

## 三、三角形与三角形相交检测：`Triangle_2_Triangle_2` do_intersect

**文件**：`Intersections_2/internal/Triangle_2_Triangle_2_do_intersect_impl.h`

**算法**：基于 Philippe Guigue 的方法，仅使用 `orientation_2` 谓词，无需计算实际交点。

**复杂度**：O(1)，最多调用约 10 次 `orientation_2`。

### 核心思想

将问题分解为：确定三角形 T1 的顶点 P1 相对于三角形 T2 的三条边的位置，然后根据 P1 的位置选择不同的子测试。

**预处理**：确保两个三角形都是逆时针方向（通过交换顶点实现，不修改原始数据）：

```cpp
if (orientation(P1, Q1, R1) != POSITIVE) { q1 = &R1; r1 = &Q1; }
if (orientation(P2, Q2, R2) != POSITIVE) { q2 = &R2; r2 = &Q2; }
```

**主逻辑**：定位 P1 相对于 T2 的三条边（p2q2、q2r2、r2p2）的位置，分为 8 种情况：

- P1 在三条边的正侧 → P1 在 T2 内部 → 相交
- P1 在某条边的负侧 → 调用 `intersection_test_edge` 或 `intersection_test_vertex`

### `intersection_test_vertex`：顶点-顶点情形

当 P1 "看到" T2 的某个顶点（即 P1 在该顶点对应的两条边的负侧）时调用。此时需要判断 T1 的边 P1Q1 或 P1R1 是否穿过 T2。

### `intersection_test_edge`：顶点-边情形

当 P1 "看到" T2 的某条边时调用。此时需要判断 T1 的边 P1Q1 是否穿过该边。

---

## 四、三角形与三角形相交计算：`Triangle_2_Triangle_2` intersection

**文件**：`Intersections_2/internal/Triangle_2_Triangle_2_intersection_impl.h`

**算法**：Sutherland-Hodgman 多边形裁剪算法。

**返回类型**：`Point_2`、`Segment_2`、`Triangle_2` 或 `std::vector<Point_2>`（凸多边形，最多 6 个顶点）。

### 核心数据结构：`Pointlist_2_`

用链表存储当前多边形的顶点，每个节点记录顶点坐标和相对于当前裁剪线的方向（`Oriented_side`）：

```cpp
struct Pointlist_2_rec_ {
  Pointlist_2_rec_ *next;
  Point_2 point;
  Oriented_side side;  // ON_POSITIVE_SIDE / ON_NEGATIVE_SIDE / ON_ORIENTED_BOUNDARY
};
```

### 算法步骤

1. **初始化**：将 T1 的三个顶点放入链表，作为初始多边形。
2. **确定 T2 的方向**：检查 T2 是否逆时针，若不是则翻转边的方向。
3. **依次用 T2 的三条边裁剪**：对每条边调用 `_cut_off`。
4. **根据剩余顶点数判断结果**：0→不相交，1→点，2→线段，3→三角形，4-6→多边形。

### `_cut_off` 函数

这是算法的核心，用一条有向直线裁剪当前多边形，保留正侧的部分：

```
对每条边 (last, cur)：
  若 last 和 cur 在直线两侧：
    计算交点，插入链表
删除所有在负侧的顶点
```

**退化处理**：当 `list_size == 2`（多边形退化为线段）且裁剪后新增了一个边界点时，需要删除该点（因为线段与直线的交点不应增加顶点数）。

### 方向修正

最终返回多边形时，检查顶点是否为顺时针顺序，若是则翻转：

```cpp
// 对精确数值类型，直接用 right_turn 判断
// 对浮点类型，用 Shoelace 公式的符号判断
if (Is_cw<K, typename Algebraic_structure_traits<typename K::FT>::Is_exact>()(points))
  std::reverse(points.begin(), points.end());
```

---

## 五、直线与直线相交：`Line_2_Line_2`

**文件**：`Intersections_2/Line_2_Line_2.h`

**算法**：直接求解线性方程组。

直线方程为 `ax + by + c = 0`。两直线 (a1,b1,c1) 和 (a2,b2,c2) 的交点满足：

```
det = a1*b2 - a2*b1
x = (b1*c2 - b2*c1) / det
y = (a2*c1 - a1*c2) / det
```

CGAL 使用齐次坐标避免除法：交点为 `Point_2(b1*c2 - b2*c1, a2*c1 - a1*c2, det)`。

**三种情况**：
- `det != 0`：唯一交点
- `det == 0, num == 0`：两直线重合，返回 `Line_2`
- `det == 0, num != 0`：两直线平行，无交点

---

## 六、矩形（Iso_rectangle_2）与线段/射线/直线相交

**文件**：`Intersections_2/Iso_rectangle_2_Segment_2.h` 等

**算法**：Cohen-Sutherland 裁剪思想的变体，通过参数化区间求交。

对于线段 P + t*(Q-P)，t ∈ [0,1]，分别计算线段与矩形四条边所在直线的交参数，求出线段在矩形内的参数区间 [t_enter, t_exit]。若区间非空，则相交。

**实现细节**：使用 `Straight_2_` 辅助类（`internal/Straight_2.h`）维护参数区间，通过 `cut_off_*` 系列函数逐步缩小区间。

---

## 七、Bbox_2 相关相交

**文件**：`Intersections_2/Bbox_2_*.h`

`Bbox_2` 使用 `double` 坐标，与 Kernel 无关。这些函数直接操作 `double` 值，不通过 Kernel 谓词，因此速度更快但不保证精确性。

`Bbox_2` 与 `Bbox_2` 的相交检测是最简单的情况：

```cpp
bool do_intersect(const Bbox_2& bb1, const Bbox_2& bb2) {
  return bb1.xmax() >= bb2.xmin() && bb2.xmax() >= bb1.xmin() &&
         bb1.ymax() >= bb2.ymin() && bb2.ymax() >= bb1.ymin();
}
```

---

## 八、`Intersection_traits_2` 的设计

`Intersection_traits_2` 通过模板特化为每对几何类型定义交集的可能类型。例如：

| 对象对 | 可能的交集类型 |
|---|---|
| `Point_2` ∩ `Point_2` | `Point_2` |
| `Segment_2` ∩ `Segment_2` | `Point_2` 或 `Segment_2` |
| `Line_2` ∩ `Line_2` | `Point_2` 或 `Line_2` |
| `Triangle_2` ∩ `Triangle_2` | `Point_2`、`Segment_2`、`Triangle_2` 或 `vector<Point_2>` |
| `Iso_rectangle_2` ∩ `Line_2` | `Point_2`、`Segment_2` 或 `Iso_rectangle_2` |

返回类型统一为 `std::optional<std::variant<T1, T2, ...>>`，通过 `intersection_return<>` 辅助函数构造。

---

## 九、关键设计模式总结

### 谓词优先，构造延迟

`do_intersect` 函数只调用谓词（`orientation_2`、`compare_xy_2` 等），不进行任何坐标计算。`intersection` 函数在确认相交后才计算交点坐标。这是 CGAL 的核心设计原则，保证了精确性。

### 端点复用

当交点恰好是某个输入端点时，直接返回该端点对象，而不是重新计算。这避免了浮点误差，也是 `S2S2_inter_info::pt_ids` 字段存在的原因。

### 静态过滤（Static Filters）

对于 `Bbox_2` 相关的相交检测，CGAL 使用静态过滤器（`use_static_filters` 模板参数）：先用浮点数快速计算，若结果不确定（在误差范围内）再用精确算术重新计算。这是 `Filtered_kernel` 的核心机制在 Intersections 包中的体现。
