# CGAL Polygon Package 算法与实现分析

## 概述

本文档分析 Polygon 包中各核心算法的实现细节，涵盖简单性检测、凸性检测、点包含测试、面积计算、方向判断等。所有算法实现位于以下文件中：

- `include/CGAL/Polygon_2_algorithms.h` — 算法声明与无 Traits 版本的包装
- `include/CGAL/Polygon_2/Polygon_2_algorithms_impl.h` — 算法实现
- `include/CGAL/Polygon_2/Polygon_2_simplicity.h` — 简单性检测的扫描线实现

---

## 一、简单性检测：`is_simple_2`

### 算法：平面扫描（Plane Sweep）

**复杂度**：O(n log n)，n 为顶点数。

**核心思想**：用一条从左到右移动的竖直扫描线，维护一棵按 y 坐标排序的活跃边集合（用 `std::set` 实现）。如果在扫描过程中发现任何非相邻边相交，则多边形不简单。

**事件类型**：扫描线在每个顶点处触发三类事件：

1. **插入事件（Insertion Event）**：顶点的两条相邻边都向右延伸。将两条边同时插入活跃边集合，并检查它们与相邻边是否相交。
2. **替换事件（Replacement Event）**：一条边向左延伸，另一条向右延伸。用右边替换左边，检查顶点是否在上下相邻边之间。
3. **删除事件（Deletion Event）**：两条相邻边都向左延伸。从集合中删除两条边，检查顶点是否在上下相邻边之间。

**关键实现细节**：

活跃边集合的比较函数 `Less_segments` 是整个算法的核心。它需要处理一个特殊情况：当一条新边（尚未在树中）与树中的边比较时，新边被视为一个"点"（其左端点），通过 `orientation_2` 谓词判断该点在树中边的哪一侧：

```cpp
// 新边的左端点 mid 与树中边 (left, right) 的方向关系
switch (orientation_2(point(left), point(mid), point(right))) {
  case LEFT_TURN: return true;   // mid 在边的左侧（下方）
  case RIGHT_TURN: return false; // mid 在边的右侧（上方）
  case COLLINEAR: break;         // 共线 → 不简单
}
```

**退化处理**：在扫描前，算法先对顶点排序并检查是否有重复顶点（度数 > 2 的情况），这是对扫描线算法的一个补丁修复：

```cpp
// 排序后检查相邻重复点
std::sort(points.begin(), points.end(), polygon_traits.less_xy_2_object());
for(; succ != points.end(); ++it, ++succ) {
  if(equal_2(*it, *succ)) return false;
}
```

---

## 二、凸性检测：`is_convex_2`

**复杂度**：O(n)

**算法**：遍历所有相邻三元组 (prev, current, next)，检查两个条件：

1. **方向一致性**：所有三元组的方向（`orientation_2`）必须一致，不能既有顺时针又有逆时针的三元组。
2. **方向变化次数**：按 x 坐标排序的相邻顶点对，方向变化（从 less 到 greater 或反之）不能超过 2 次。超过 2 次说明多边形不是凸的。

```cpp
bool HasClockwiseTriples = false;
bool HasCounterClockwiseTriples = false;
int NumOrderChanges = 0;

// 遍历所有三元组
switch (orientation(*previous, *current, *next)) {
  case CLOCKWISE:        HasClockwiseTriples = true; break;
  case COUNTERCLOCKWISE: HasCounterClockwiseTriples = true; break;
  case ZERO: /* 处理共线点 */ break;
}

if (HasClockwiseTriples && HasCounterClockwiseTriples) return false;
if (NumOrderChanges > 2) return false;
```

**共线点处理**：算法通过 `Equal_2` 谓词跳过重复点，通过 `goto` 语句在遇到共线三元组时重新检查下一个三元组（这是代码中少见的 `goto` 用法，用于处理退化情况）。

---

## 三、点包含测试：`bounded_side_2`

**复杂度**：O(n)

**算法**：射线投射法（Ray Casting）。从查询点向右发射一条水平射线，统计与多边形边的交叉次数。奇数次交叉 → 点在内部；偶数次 → 点在外部。

**关键退化处理**：射线恰好穿过顶点时的处理是实现中最复杂的部分。CGAL 采用"半开区间"约定：

- 对于非水平边，**上端点属于该边**，下端点不属于。
- 这保证了每个顶点只被计数一次，避免了顶点被重复计数导致奇偶性错误。

实现通过一个嵌套的 switch 语句处理所有情况（当前顶点 y 与查询点 y 的比较结果 × 下一顶点 y 与查询点 y 的比较结果）：

```cpp
// 辅助函数：判断点在线段 (low, high) 的哪一侧
// 先用 x 坐标快速过滤，再用 orientation_2 精确判断
int which_side_in_slab(point, low, high, orientation_2, compare_x_2)
```

`which_side_in_slab` 函数先尝试用 x 坐标比较快速判断（避免调用昂贵的 `orientation_2`），只有在 x 坐标无法区分时才调用 `orientation_2`，这是一个性能优化。

---

## 四、方向判断：`orientation_2`

**复杂度**：O(n)（需要先找最左顶点）

**算法**：找到多边形的最左顶点（x 坐标最小，平局取 y 最小），然后对该顶点及其前后相邻顶点组成的三元组调用 `orientation_2` 谓词。

**正确性依据**：最左顶点处的局部方向等于整个简单多边形的全局方向。这是因为最左顶点是一个"极值顶点"，其局部凸性与全局方向一致。

**前置条件**：要求多边形是简单的（`CGAL_precondition(is_simple_2(...))`）。内部还提供了一个无前置条件版本 `orientation_2_no_precondition`，供 Straight Skeleton 等包在处理非严格简单多边形时使用。

---

## 五、面积计算：`polygon_area_2` / `area_2`

**复杂度**：O(n)

**算法**：Shoelace 公式（高斯面积公式）的三角剖分变体。以第一个顶点为基点，将多边形分解为 n-2 个三角形，累加各三角形的有向面积：

```cpp
ForwardIterator third = second;
while (++third != last) {
    result = result + compute_area_2(*first, *second, *third);
    second = third;
}
```

`compute_area_2` 计算三点组成的有向三角形面积（即行列式的一半）。对于凸多边形，所有三角形面积同号；对于非凸简单多边形，部分三角形面积为负，但总和仍然正确。

**注意**：对于非简单多边形，面积无明确定义，但函数仍会返回一个值（自相交区域会被抵消）。

---

## 六、极值顶点查找：`left/right/top/bottom_vertex_2`

**复杂度**：O(n)

这四个函数都是对 `std::min_element` / `std::max_element` 的简单包装，使用 Traits 提供的比较函数对象：

```cpp
// left_vertex_2: 最小 x，平局取最小 y
return std::min_element(first, last, traits.less_xy_2_object());

// top_vertex_2: 最大 y，平局取最大 x
return std::max_element(first, last, traits.less_yx_2_object());
```

`Less_xy_2` 按 (x, y) 字典序比较，`Less_yx_2` 按 (y, x) 字典序比较，这两个比较器由 Kernel 提供。

---

## 七、有向侧判断：`oriented_side_2`

**复杂度**：O(n)

这是 `bounded_side_2` 和 `orientation_2` 的组合：

```cpp
Orientation o = orientation_2(first, last, traits);  // 多边形方向
Bounded_side b = bounded_side_2(first, last, point, traits);  // 点的位置

// 根据多边形方向将 bounded_side 转换为 oriented_side
if (b == ON_BOUNDED_SIDE)
    return (o == CLOCKWISE) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
```

对于顺时针多边形，内部是负侧；对于逆时针多边形，内部是正侧。

---

## 八、共线点过滤：`filter_collinear_points`（内部工具函数）

这是一个内部工具函数（位于 `internal::Polygon_2` 命名空间），用于在简单性检测前预处理多边形，去除近似共线的连续顶点。

**算法**：滑动窗口，对每个连续三元组 (o, p, q) 计算行列式：

```
det = |o.x - q.x  o.y - q.y|
      |p.x - q.x  p.y - q.y|
```

如果 `|det| <= tolerance`，则跳过中间点 p（认为三点共线）。

---

## 九、迭代器与循环子的实现

### `Polygon_2_edge_iterator`

边迭代器不存储边，而是存储顶点迭代器，按需构造 `Segment_2`：

```cpp
// 解引用时动态构造边
Segment_2 operator*() const {
    return traits->construct_segment_2_object()(*m_it, *next_it);
}
```

### `Polygon_2_vertex_circulator`

顶点循环子基于 CGAL 的 `Circulator` 机制，在到达容器末尾时自动回绕到开头，实现循环遍历语义。

---

## 十、数据结构层次的实现

### `General_polygon_with_holes_2<Polygon_>`

- 外边界：单个 `Polygon_` 对象（`m_pgn`）
- 洞：`std::deque<Polygon_>` 容器（`m_holes`）

使用 `deque` 而非 `vector` 是因为洞的数量通常较少，且 `deque` 在两端插入删除更高效。

### `Multipolygon_with_holes_2`

- 内部存储：`std::deque<Polygon_with_holes_2>`

提供了 `begin()`/`end()` 接口，使其可以直接用于范围 for 循环。

---

## 十一、等价性判断：`operator==`

`Polygon_2` 的相等判断考虑了**循环等价**：两个多边形相等当且仅当存在一个循环置换使得顶点序列相同。实现在 `Polygon_2_impl.h` 中，通过在第二个多边形中搜索第一个多边形的起始顶点，然后逐一比较。

`Polygon_with_holes_2` 的相等判断：外边界相等，且洞的集合相等（不考虑洞的顺序，通过在临时列表中查找并删除实现）。
