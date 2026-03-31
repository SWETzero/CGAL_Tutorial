# CGAL Polygon Package 设计哲学

## 概述

Polygon 包是 CGAL 中最基础的几何包之一，提供二维多边形的表示与操作。它的设计体现了 CGAL 整体的核心哲学：**泛型编程、算法与数据结构分离、精确性与效率的平衡**。

---

## 一、泛型编程：Kernel 参数化

`Polygon_2` 是一个模板类，接受两个模板参数：

```cpp
template <class Traits_, class Container_ = std::vector<typename Traits_::Point_2>>
class Polygon_2;
```

- `Traits_`（通常是一个 Kernel）决定了点的类型、坐标类型（`FT`）以及所有几何谓词（方向判断、面积计算等）的具体实现。
- `Container_` 决定了顶点的存储方式，默认为 `std::vector`，但可以替换为任何满足 STL 容器要求的类型。

这种设计的核心思想是：**几何算法的正确性由 Kernel 保证，而不是硬编码在算法里**。用户可以选择 `Exact_predicates_inexact_constructions_kernel` 获得精确谓词，或者选择其他 Kernel 在精度和性能之间权衡，而多边形的所有操作代码无需改变。

---

## 二、算法与数据结构分离

Polygon 包将**算法**和**数据结构**明确分开：

- `Polygon_2_algorithms.h` 提供了一组作用于**点序列迭代器**的全局函数（`is_simple_2`、`bounded_side_2`、`orientation_2` 等）。
- `Polygon_2` 类是对这些全局函数的封装，它的成员函数（如 `is_simple()`、`bounded_side()`）本质上只是将内部容器的迭代器传递给对应的全局函数。

这意味着：**你不需要使用 `Polygon_2` 类才能使用这些算法**。任何持有点序列的代码都可以直接调用全局函数。这是 CGAL 中"算法作用于迭代器范围"这一设计原则的典型体现，与 STL 的设计哲学一脉相承。

---

## 三、无缓存原则（No Caching）

`Polygon_2` 是一个轻量级的容器包装器，**不缓存任何计算结果**。每次调用 `is_simple()`、`area()`、`orientation()` 等方法，都会重新计算。

文档中明确说明了这一点：

> The `Polygon_2` class is a wrapper around a container of points, but little more. Especially, computed values are not cached.

这个设计决策的背后逻辑是：
- 多边形的顶点可以随时被修改（插入、删除、移动），缓存失效的管理会引入复杂性。
- 将缓存策略的决定权交给用户，而不是在库内部强制一种策略。
- 保持类的职责单一：`Polygon_2` 只负责存储和提供访问接口，不负责维护派生属性。

这是一种**简单性优先**的设计哲学，代价是重复调用同一谓词时有额外开销，但避免了状态管理的复杂性。

---

## 四、概念（Concept）驱动的接口设计

Polygon 包定义了若干 Concept：

- `PolygonTraits_2`：规定了 Traits 类必须提供的类型和函数对象（`Less_xy_2`、`Orientation_2`、`Compute_area_2` 等）。
- `GeneralPolygonWithHoles_2`：规定了带洞多边形类必须满足的接口。
- `MultipolygonWithHoles_2`：规定了多多边形类的接口。

这种 Concept 驱动的设计使得：
- 算法只依赖 Concept，而不依赖具体类型，实现了真正的泛型。
- 用户可以提供自定义的 Traits 类，只要满足 `PolygonTraits_2` 的要求，所有算法都能正常工作。
- 例如，`Projection_traits_xy_3` 可以作为 Traits 传入，让 2D 多边形算法直接作用于 3D 空间中的投影多边形。

---

## 五、层次化的数据模型

Polygon 包提供了三个层次的数据结构，形成清晰的组合关系：

```
Multipolygon_with_holes_2
    └── Polygon_with_holes_2  (多个)
            ├── outer_boundary: Polygon_2
            └── holes: Polygon_2  (多个)
```

- `General_polygon_with_holes_2<Polygon_>` 是泛型基类，接受任意多边形类型作为模板参数，用于支持非线性多边形（如曲线多边形）的场景（Boolean Set Operations 包就使用了这个泛型版本）。
- `Polygon_with_holes_2<Kernel, Container_>` 是具体化版本，将 `Polygon_2` 作为边界类型。
- `Multipolygon_with_holes_2` 是最顶层的容器，持有多个 `Polygon_with_holes_2`。

这种设计体现了**组合优于继承**的原则，同时通过泛型基类支持了更广泛的使用场景。

---

## 六、STL 兼容性

`Polygon_2` 的设计高度兼容 STL：

- 顶点迭代器（`Vertex_iterator`）和边迭代器（`Edge_const_iterator`）都满足 STL 迭代器要求。
- 同时提供了 CGAL 特有的**循环子（Circulator）**接口（`Vertex_circulator`、`Edge_const_circulator`），用于表达多边形的循环拓扑结构。
- 提供了 `vertices()` 和 `edges()` 方法返回范围（`Iterator_range`），支持现代 C++ 的范围 for 循环。
- 支持 `begin()`/`end()` 接口，使 `Polygon_2` 本身可以作为范围使用。

这种设计让 `Polygon_2` 可以无缝接入 STL 算法生态，同时通过 Circulator 表达了多边形特有的循环结构语义。

---

## 七、精确性保证的责任分配

Polygon 包的算法（如 `is_simple_2`、`bounded_side_2`）在实现上对数值退化情况（如射线恰好穿过顶点、三点共线等）都有明确的处理逻辑，但**精确性的最终保证由 Kernel 提供**。

文档中对 `is_simple_2` 的说明：

> The algorithm is quite robust when used with inexact number types.

这表明算法本身在设计上考虑了鲁棒性，但如果需要绝对精确的结果，应使用精确 Kernel（如 `Exact_predicates_exact_constructions_kernel`）。这种责任分配方式是 CGAL 的一贯做法：算法尽量鲁棒，精确性由 Kernel 层保证。
