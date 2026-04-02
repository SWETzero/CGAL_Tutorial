# CGAL Intersections_3 Package 算法与实现分析

## 概述

本文档分析 Intersections_3 包中各核心算法的实现细节，涵盖 3D 几何对象之间的相交检测（`do_intersect`）与相交计算（`intersection`）。所有算法实现位于以下文件中：

- `include/CGAL/intersection_3.h` — 公开 API 入口
- `include/CGAL/Intersection_traits_3.h` — 返回类型 traits
- `include/CGAL/Intersections_3/XxxType_YyyType.h` — 各对象对的薄包装层
- `include/CGAL/Intersections_3/internal/` — 所有实际算法实现

---

## 一、架构设计

Intersections_3 的架构与 Intersections_2 完全一致，但 3D 情况下几何对象更多（增加了 `Plane_3`、`Sphere_3`、`Tetrahedron_3`、`Iso_cuboid_3` 等），对象对的数量也更多。

每个 `XxxType_YyyType.h` 文件是一个薄包装层，仅包含 `#include` 和一个转发调用：

```cpp
// Plane_3_Triangle_3.h
#include <CGAL/Intersections_3/internal/Plane_3_Triangle_3_do_intersect.h>
#include <CGAL/Intersections_3/internal/Plane_3_Triangle_3_intersection.h>
// ...
CGAL_INTERSECTION_FUNCTION(Plane_3, Triangle_3, 3)
CGAL_DO_INTERSECT_FUNCTION(Plane_3, Triangle_3, 3)
```

实际算法全部在 `internal/` 子目录中，文件名格式为 `TypeA_TypeB_do_intersect.h` 和 `TypeA_TypeB_intersection.h`。

---

## 二、直线与平面相交：`Line_3_Plane_3`

**文件**：`internal/Line_3_Plane_3_intersection.h`

**返回类型**：`Point_3` 或 `Line_3`（直线在平面内时）

### 算法

平面方程：`ax + by + cz + d = 0`。直线用点 `P` 和方向 `D` 表示：`P + t*D`。

代入平面方程求参数 t：

```
t = -(a*Px + b*Py + c*Pz + d) / (a*Dx + b*Dy + c*Dz)
```

CGAL 使用齐次坐标（`hx, hy, hz, hw`）避免除法，直接构造交点：

```cpp
RT num = plane.a()*line_pt.hx() + plane.b()*line_pt.hy()
       + plane.c()*line_pt.hz() + wmult_hw(plane.d(), line_pt);
RT den = plane.a()*line_dir.dx() + plane.b()*line_dir.dy()
       + plane.c()*line_dir.dz();

// 交点（齐次坐标）：
Point_3(den*line_pt.hx() - num*line_dir.dx(),
        den*line_pt.hy() - num*line_dir.dy(),
        den*line_pt.hz() - num*line_dir.dz(),
        wmult_hw(den, line_pt))
```

**三种情况**：
- `den != 0`：唯一交点
- `den == 0, num == 0`：直线在平面内，返回 `Line_3`
- `den == 0, num != 0`：直线平行于平面，无交点

`wmult_hw` 是一个辅助函数，处理齐次坐标的 `hw` 分量乘法，对 Cartesian kernel（`hw=1`）直接返回参数，对 Homogeneous kernel 则乘以 `hw`。

---

## 三、三角形与三角形相交检测：`Triangle_3_Triangle_3` do_intersect

**文件**：`internal/Triangle_3_Triangle_3_do_intersect.h`

**算法**：基于 Philippe Guigue 和 Olivier Devillers 的论文《A Fast and Robust Triangle-Triangle Intersection Test》，仅使用 `orientation_3` 和 `coplanar_orientation_3` 谓词。

**复杂度**：O(1)，最多约 15 次谓词调用。

### 核心思想：区间重叠测试

两个三角形 T1(p,q,r) 和 T2(a,b,c) 相交，当且仅当它们的交线段在两个三角形的平面投影上重叠。

**步骤一**：计算 T1 的三个顶点相对于 T2 所在平面的方向：

```cpp
Orientation dp = orientation(a, b, c, p);  // p 在 T2 平面的哪侧
Orientation dq = orientation(a, b, c, q);
Orientation dr = orientation(a, b, c, r);
```

若三个顶点同侧（全正或全负），则不相交，直接返回 false。

**步骤二**：确定"孤立顶点"。T1 与 T2 平面的交线将 T1 的三条边分为两组：两条边与平面相交（形成区间端点），一条边不相交。孤立顶点是与另外两个顶点在平面不同侧的那个顶点。

例如，若 dp=POSITIVE，dq=NEGATIVE，dr=NEGATIVE，则 p 是孤立顶点，交线段的端点在边 pq 和 rp 上：

```cpp
s_min1 = &p; t_min1 = &q;  // 边 pq 上的交点
s_max1 = &r; t_max1 = &p;  // 边 rp 上的交点
```

**步骤三**：对 T2 重复上述过程，得到 T2 与 T1 平面的交线段端点 (s_min2, t_min2) 和 (s_max2, t_max2)。

**步骤四**：检查两个区间是否重叠：

```cpp
return (orientation(s_min1, t_min1, s_min2, t_min2) != POSITIVE) &&
       (orientation(s_max1, t_max1, t_max2, s_max2) != POSITIVE);
```

这两个 `orientation_3` 调用等价于检查两个有向区间是否有重叠部分，无需计算实际的交线坐标。

### 共面情形：`do_intersect_coplanar`

当 T1 的三个顶点都在 T2 的平面上时（dp=dq=dr=COPLANAR），退化为 2D 问题。

算法与 `Triangle_2_Triangle_2` 的 `do_intersect` 完全类似，但使用 `coplanar_orientation_3`（在已知共面的情况下，用 3D 点计算 2D 方向）代替 `orientation_2`。

定位 T1 的顶点 p 相对于 T2 的三条边的位置，然后调用 `_intersection_test_edge` 或 `_intersection_test_vertex`。

---

## 四、三角形与三角形相交计算：`Triangle_3_Triangle_3` intersection

**文件**：`internal/Triangle_3_Triangle_3_intersection.h`

**返回类型**：`Point_3`、`Segment_3`、`Triangle_3` 或 `std::vector<Point_3>`（共面时最多 6 个顶点的多边形）

### 非共面情形

**算法**：

1. 计算两平面的交线 L（调用 `Plane_3_Plane_3_intersection`）。
2. 计算 T1 与 L 的交集（调用 `Line_3_Triangle_3_intersection`），得到线段或点 I1。
3. 计算 T2 与 L 的交集，得到 I2。
4. 计算 I1 和 I2 的交集（两个共线线段的交集）。

### 共面情形

退化为 2D 问题，使用 Sutherland-Hodgman 裁剪（与 `Triangle_2_Triangle_2` 的 intersection 相同思路）。

### `Point_on_triangle` 结构体

这是实现中最复杂的部分。为了精确处理退化情形（顶点在对方三角形的边上等），每个交点用 `Point_on_triangle` 结构体表示，记录该点在两个三角形上的"身份"：

```cpp
struct Point_on_triangle {
  std::pair<int, int> t1_t2_ids;
  // t1_t2_ids.first: 点在 T1 上的位置
  //   -1,-2,-3: T1 的顶点 p,q,r
  //   1,2,3: T1 的边 pq,qr,rp 上的点
  //   0: T2 的顶点在 T1 内部
  // t1_t2_ids.second: 类似地表示在 T2 上的位置
  FT alpha;  // 若在边上，alpha 表示在边上的参数位置
};
```

这个设计使得在判断点的方向关系时，可以直接用组合信息推断，而不需要重新计算 `orientation`，避免了精度问题。

---

## 五、Bbox_3 与线段相交检测：`Bbox_3_Segment_3`

**文件**：`internal/Bbox_3_Segment_3_do_intersect.h`

**算法**：参数化区间法（Amy Williams 等人的"slab method"）。

**核心思想**：将线段参数化为 P + t*(Q-P)，t ∈ [0,1]。分别计算线段与 x、y、z 三个轴对齐平板（slab）的参数区间，求三个区间的交集。若交集与 [0,1] 有重叠，则相交。

### 实现细节

对每个坐标轴（以 x 为例），计算线段进入和离开 bbox 的参数值：

```cpp
// 若 qx >= px（线段向右）
tmin = (bxmin - px) / (qx - px)  // 进入 x-slab 的参数
tmax = (bxmax - px) / (qx - px)  // 离开 x-slab 的参数
```

为避免除法，CGAL 用分子/分母的形式维护区间端点，通过交叉乘法比较：

```
t1/d1 > t2/d2  ⟺  t1*d2 > t2*d1  （当 d1,d2 > 0 时）
```

### 静态过滤器：`Do_intersect_bbox_segment_aux_is_greater`

这是一个精心设计的过滤器，有两个特化版本：

**通用版本**（非 double）：直接比较，无过滤。

**double 特化版本**（`use_static_filters=true`）：
1. 在处理每个坐标轴时，记录输入值的最大绝对值（`tmax`, `dmax`）。
2. 计算误差界：`error = tmax * dmax * EPS`（EPS ≈ 8.89e-16，约为 4 个 ulp）。
3. 比较 `a - b` 与误差界：若 `|a-b| > error`，结果确定；否则返回 `Uncertain::indeterminate()`，触发精确计算。

```cpp
result_type operator()(const FT& a, const FT& b) const {
  const FT x = a - b;
  if (x > error)  return true;
  if (x < -error) return false;
  return uncertain();  // 触发精确重算
}
```

**溢出/下溢保护**：若 `dmax > 1e153` 或 `dmax < 1e-146`，跳过静态过滤直接用精确算术。

### 模板参数 `bounded_0` 和 `bounded_1`

- `bounded_0=true`：线段的起点有界（即不是射线/直线的起点方向无穷远）
- `bounded_1=true`：线段的终点有界

这使得同一套代码可以处理线段（bounded_0=true, bounded_1=true）、射线（bounded_0=true, bounded_1=false）和直线（bounded_0=false, bounded_1=false）。

---

## 六、Bbox_3 与三角形相交检测：`Bbox_3_Triangle_3`

**文件**：`internal/Bbox_3_Triangle_3_do_intersect.h`

**算法**：Separating Axis Theorem（SAT，分离轴定理）。

**核心思想**：两个凸体不相交，当且仅当存在一个分离轴，使得两个凸体在该轴上的投影不重叠。对于 Bbox 和三角形，需要测试以下轴：

1. **3 个坐标轴**（x, y, z）：检查三角形的 bbox 与 Bbox 是否重叠。
2. **三角形法向量**：检查 Bbox 的 8 个顶点是否都在三角形平面的同侧。
3. **9 个边叉积轴**：三角形的 3 条边 × Bbox 的 3 个坐标轴方向 = 9 个轴。

**实现**：代码展开了所有 9 个边叉积测试，每个测试形如：

```cpp
// 测试轴 e0 × ex（e0 是三角形第一条边，ex 是 x 轴方向）
// 叉积 = (0, -e0z, e0y)
// 三角形顶点在该轴上的投影：
FT p0 = -e0.z()*v0.y() + e0.y()*v0.z();
FT p1 = -e0.z()*v1.y() + e0.y()*v1.z();
FT p2 = -e0.z()*v2.y() + e0.y()*v2.z();
// Bbox 在该轴上的半径：
FT r = half_x * CGAL::abs(e0.z()) + half_y * CGAL::abs(e0.y());
// 分离条件：
if (CGAL::min(p0,p1,p2) > r || CGAL::max(p0,p1,p2) < -r) return false;
```

---

## 七、平面与平面相交：`Plane_3_Plane_3`

**文件**：`internal/Plane_3_Plane_3_intersection.h`

**返回类型**：`Line_3` 或 `Plane_3`（两平面重合时）

### 算法

两平面 (a1,b1,c1,d1) 和 (a2,b2,c2,d2) 的交线方向为法向量的叉积：

```
direction = (b1*c2 - b2*c1, a2*c1 - a1*c2, a1*b2 - a2*b1)
```

交线上的一个点通过求解线性方程组得到。CGAL 选择叉积中绝对值最大的分量对应的坐标轴，将该坐标设为 0，求解另外两个坐标：

```cpp
// 若 |b1*c2 - b2*c1| 最大，令 x=0，求解 y 和 z
// b1*y + c1*z = -d1
// b2*y + c2*z = -d2
```

这种选择保证了数值稳定性（避免用接近零的分量做除数）。

---

## 八、三平面相交：`Plane_3_Plane_3_Plane_3`

**文件**：`internal/Plane_3_Plane_3_Plane_3_intersection.h`

**返回类型**：`Point_3`、`Line_3` 或 `Plane_3`

**算法**：Cramer 法则求解 3×3 线性方程组：

```
a1*x + b1*y + c1*z = -d1
a2*x + b2*y + c2*z = -d2
a3*x + b3*y + c3*z = -d3
```

行列式 `det = a1*(b2*c3 - b3*c2) - b1*(a2*c3 - a3*c2) + c1*(a2*b3 - a3*b2)`。

- `det != 0`：唯一交点（用 Cramer 法则计算）
- `det == 0`：退化情形，需进一步判断是交线还是无交集

---

## 九、四面体相关相交

**文件**：`internal/Tetrahedron_3_Bounded_3_do_intersect.h`、`internal/Tetrahedron_3_Unbounded_3_do_intersect.h`

四面体（Tetrahedron_3）的相交检测通过将其分解为 4 个三角面来实现。

**`Tetrahedron_3_Bounded_3_do_intersect`**（有界对象：点、线段、三角形、另一个四面体）：

对于点：检查点是否在四面体的 4 个半空间内（即在每个面的正侧）。

对于线段/三角形：检查是否与四面体的任意一个面相交，或者对象是否完全在四面体内部。

**`Tetrahedron_3_Unbounded_3_do_intersect`**（无界对象：直线、射线、平面）：

对于直线/射线：计算与四面体的 4 个面的交参数，求参数区间的交集。

---

## 十、球体相关相交

### `Sphere_3_Sphere_3`

**文件**：`internal/Sphere_3_Sphere_3_intersection.h`

**返回类型**：`Point_3`（外切）、`Circle_3`（相交圆）或 `Sphere_3`（重合）

**算法**：两球心距离 d 与半径之和/差的比较：
- `d > r1 + r2`：不相交
- `d == r1 + r2`：外切，交点在连心线上
- `|r1 - r2| < d < r1 + r2`：相交圆，交圆所在平面垂直于连心线
- `d == |r1 - r2|`：内切
- `d == 0, r1 == r2`：重合

交圆的圆心和半径通过几何关系计算：

```
// 交圆圆心到球1圆心的距离
h = (d² + r1² - r2²) / (2*d)
// 交圆半径
r = sqrt(r1² - h²)
```

### `Line_3_Sphere_3` do_intersect

**文件**：`internal/Line_3_Sphere_3_do_intersect.h`

**算法**：计算直线到球心的距离的平方，与半径平方比较：

```
dist² = ||(P - C) - ((P-C)·D)*D||²
```

其中 P 是直线上的点，C 是球心，D 是直线方向（单位向量）。

CGAL 用精确算术计算，避免开平方：

```cpp
// 判别式 Δ = r² - dist²
// Δ > 0: 两个交点
// Δ == 0: 一个切点
// Δ < 0: 不相交
```

---

## 十一、Iso_cuboid_3（轴对齐包围盒）相关相交

`Iso_cuboid_3` 是精确的轴对齐长方体（与 `Bbox_3` 不同，它使用 Kernel 的数值类型）。

### `Iso_cuboid_3_Segment_3` intersection

**文件**：`internal/Iso_cuboid_3_Segment_3_intersection.h`

**算法**：与 `Bbox_3_Segment_3` 类似的参数化区间法，但使用 Kernel 谓词而非 double 比较，保证精确性。

### `Iso_cuboid_3_Plane_3` intersection

**文件**：`internal/Iso_cuboid_3_Plane_3_intersection.h`

**算法**：计算长方体的 8 个顶点相对于平面的方向，收集在平面正侧和负侧的顶点，然后对每条穿越平面的边计算交点，构成交多边形（最多 6 个顶点的凸多边形）。

---

## 十二、关键设计模式总结

### do_intersect 与 intersection 的分离

所有对象对都严格分离 `do_intersect`（仅用谓词）和 `intersection`（计算坐标）。这使得在只需要判断是否相交时（如碰撞检测），可以避免昂贵的坐标计算。

### 齐次坐标的使用

3D 交点计算大量使用齐次坐标（`hx, hy, hz, hw`），避免了中间除法，保证了精确性。对于 Cartesian kernel，`hw=1`，这些计算退化为普通坐标运算。

### 共面退化的处理

3D 算法中最复杂的部分是处理共面退化情形（如两个三角形共面）。CGAL 通过检测 `orientation_3` 返回 `COPLANAR` 来识别这种情况，然后切换到 2D 算法（使用 `coplanar_orientation_3`）。

### 模板化的 Bbox 算法

`Bbox_3_Segment_3_do_intersect.h` 中的 `do_intersect_bbox_segment_aux` 函数通过模板参数 `bounded_0`、`bounded_1` 和 `use_static_filters` 统一处理线段、射线、直线与 Bbox 的相交，以及是否使用静态过滤器，避免了代码重复。
