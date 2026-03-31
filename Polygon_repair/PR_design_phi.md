# CGAL Polygon Repair Package 设计哲学

## 概述

Polygon Repair 包提供二维多边形的修复功能：给定可能无效的多边形输入（自相交、洞嵌套错误、多边形重叠等），输出一个有效的 `Multipolygon_with_holes_2`。它的设计体现了 CGAL 的核心哲学，同时针对"修复"这一特殊语义引入了独特的设计决策。

---

## 一、以三角剖分为核心的修复范式

Polygon Repair 的核心思想与 Polygon 包的算法（如 `is_simple_2`、`bounded_side_2`）截然不同。Polygon 包的算法假设输入是有效的，在有效输入上计算属性；而 Polygon Repair 的出发点是：**输入可能是任意无效的，修复算法必须对所有输入都能产生有效输出**。

实现这一目标的核心手段是**约束 Delaunay 三角剖分（CDT）**：

1. 将输入多边形的所有边作为约束插入 CDT。
2. 对三角剖分的每个面打标签（内部/外部/洞）。
3. 从标签重建有效的多边形。

这种"先三角化，再重建"的范式天然地处理了所有几何退化情况：自相交的边在三角化时被自动分割，重叠区域在标签阶段被正确处理，最终重建的多边形保证有效。

---

## 二、标签规则作为策略模式（Strategy Pattern）

修复算法的第二步（打标签）决定了哪些三角形属于"内部"，这直接决定了修复结果的语义。Polygon Repair 将标签规则设计为可替换的策略，通过 Tag 类传入：

```cpp
// 四种规则，均为空的 Tag 类
struct Even_odd_rule {};
struct Non_zero_rule {};
struct Union_rule {};
struct Intersection_rule {};

// 通过模板参数选择规则
Multipolygon_with_holes_2 result = CGAL::Polygon_repair::repair(polygon, Even_odd_rule());
Multipolygon_with_holes_2 result = CGAL::Polygon_repair::repair(polygon, Non_zero_rule());
```

这四种规则对应不同的修复语义：

- **Even-odd（奇偶规则）**：穿越每条边时切换内外状态，不区分边的方向。适合修复一般无效多边形。
- **Non-zero（非零绕数规则）**：根据绕数（winding number）判断内外，考虑边的方向。对同向嵌套多边形与奇偶规则结果不同。
- **Union（并集规则）**：结果为所有输入多边形的并集，适合合并多个有效多边形。
- **Intersection（交集规则）**：结果为所有输入多边形的交集。

这种设计的哲学是：**修复的"正确"结果取决于用户对输入数据的语义理解**，库不强制一种解释，而是提供多种规则供用户选择。

---

## 三、三种后端实现的分工

四种规则在实现上分为三个后端，这是 Polygon Repair 内部架构的核心：

**后端一：`Polygon_repair` 类（Even-odd 规则）**

使用 `Triangulation_with_even_odd_constraints_2` 包装的 CDT。插入约束时采用特殊的 `even_odd_insert_constraint`：如果一条边已经存在于三角化中，则**删除**它（而不是重复插入），实现奇偶抵消语义。标签阶段通过 BFS 从无穷面出发，交替标记内外区域。

**后端二：`Winding` 类（Non-zero 规则）**

使用 `Constrained_triangulation_plus_2`（CDT+），它能追踪每条约束边的方向和重数。标签阶段计算每个三角形的绕数（winding number）：从无穷面（绕数=0）出发，穿越每条约束边时根据边的方向增减绕数（`delta = +1` 或 `-1`）。绕数非零的三角形属于内部。

**后端三：`Boolean` 类（Union/Intersection 规则）**

同样使用 CDT+，但标签存储的是"覆盖层数"（`layers`）：每个三角形被多少个输入多边形覆盖。Union 规则选取 `layers > 0` 的三角形，Intersection 规则选取 `layers == N`（N 为输入多边形数量）的三角形。这使得 Union/Intersection 规则天然支持任意数量的输入多边形。

---

## 四、输出的确定性规范化

Polygon Repair 对输出施加了比"有效性"更严格的规范化约束，以保证**相同输入总产生相同输出**（确定性）：

- 相邻共线边合并（消除度数为 2 的共线顶点）
- 每条边界从字典序最小顶点开始
- 外边界逆时针，内边界（洞）顺时针
- 洞按字典序排列
- 多边形按字典序排列

这些规范化在 `reconstruct_ring` 和 `reconstruct_multipolygon` 中实现。`reconstruct_ring` 在重建环时通过线性扫描找到字典序最小顶点并旋转到首位；`reconstruct_multipolygon` 使用 `Polygon_less` 和 `Polygon_with_holes_less` 比较器将结果存入有序集合（`std::set`）。

这与 Polygon 包的"无缓存原则"形成对比：Polygon 包不对存储顺序做任何规范化，而 Polygon Repair 为了确定性输出主动施加了额外的排序开销。

---

## 五、`is_valid` 作为修复的逆操作

Polygon Repair 包提供了 `is_valid` 函数，用于验证多边形是否满足有效性定义。`is_valid` 的实现本身也使用了 CDT，这体现了一种设计上的对称性：

- `repair` 用 CDT 将无效输入转化为有效输出。
- `is_valid` 用 CDT 验证输入是否满足有效性条件。

`is_valid` 对 `Polygon_with_holes_2` 的验证流程：
1. 验证外边界和每个洞各自是简单多边形。
2. 将外边界插入 CDT，检查是否有自相交（通过 `No_constraint_intersection_tag` 捕获异常）。
3. 将洞插入 CDT，检查洞是否完全在外边界内部（通过标签传播验证）。
4. 检查内部是否连通（再次标签传播，验证只有一个内部区域）。

这种用 CDT 做验证的方式，与用 CDT 做修复的方式共享了大量基础设施（`label_region` 函数），体现了代码复用的设计意图。

---

## 六、与 Polygon 包的设计对比

| 维度 | Polygon 包 | Polygon Repair 包 |
|---|---|---|
| 输入假设 | 输入有效 | 输入可能无效 |
| 核心数据结构 | 顶点容器 + 迭代器 | 约束 Delaunay 三角剖分 |
| 算法策略 | 固定算法 | 可替换规则（Tag 类） |
| 输出类型 | 与输入类型相同 | 始终为 `Multipolygon_with_holes_2` |
| 输出规范化 | 无 | 强制字典序规范化 |
| 精确性机制 | Kernel 参数化 | Kernel 参数化 + 浮点/精确自动选择 |

两个包都依赖 Kernel 参数化，但 Polygon Repair 在 Kernel 选择上有一个额外的自动化机制：根据 `Kernel::FT` 是否为浮点类型，自动选择 `Exact_predicates_tag`（浮点时）或 `Exact_intersections_tag`（精确数值时）作为 CDT 的相交处理策略。这使得用户使用 EPICK 时也能获得正确的相交处理，无需手动指定。
