# CGAL Polygon Mesh Processing Package 设计哲学

## 概述

Polygon Mesh Processing（PMP）包是 CGAL 中最大的算法包之一，提供从基础度量计算到高级网格修复、重网格化的全套三维多边形网格处理算法。它的设计体现了 CGAL 的核心哲学，同时针对三维网格处理的特殊需求引入了若干独特的设计决策。

---

## 一、以 BGL 概念为核心的泛型设计

PMP 包的所有函数都不依赖具体的网格类型，而是依赖 Boost Graph Library（BGL）定义的图概念：

- `HalfedgeGraph`：提供半边遍历能力，用于边长计算、法向量计算等。
- `FaceGraph`：在 `HalfedgeGraph` 基础上增加面的概念，用于面积、体积计算。
- `MutableFaceGraph`：支持拓扑修改，用于重网格化、修复等操作。

这意味着 `CGAL::Surface_mesh`、`CGAL::Polyhedron_3`、OpenMesh 的 `TriMesh` 等都可以直接传入 PMP 的函数，无需任何适配代码。这是 PMP 与 Polygon 包最大的设计差异：Polygon 包依赖 Traits 参数化，而 PMP 依赖 BGL 图概念参数化。

```cpp
// 同一个函数，接受不同网格类型
CGAL::Surface_mesh<Point_3> sm;
CGAL::Polyhedron_3<K> poly;

PMP::compute_face_normals(sm, get(CGAL::face_normal, sm));   // Surface_mesh
PMP::compute_face_normals(poly, face_normals_map);            // Polyhedron_3
```

---

## 二、Named Parameters：可选参数的优雅处理

PMP 大量使用 Named Parameters 模式处理可选参数，这是与 Polygon 包最显著的 API 差异。

Polygon 包的函数签名通常是固定的（`is_simple_2(first, last, traits)`），而 PMP 的函数往往有十几个可选参数：

```cpp
// 不传任何可选参数，使用默认值
PMP::isotropic_remeshing(faces(mesh), target_length, mesh);

// 传入部分可选参数
PMP::isotropic_remeshing(faces(mesh), target_length, mesh,
    CGAL::parameters::number_of_iterations(5)
                     .protect_constraints(true)
                     .vertex_point_map(vpm));
```

Named Parameters 的实现基于 `CGAL::Named_function_parameters` 模板链，每个参数通过 `.param_name(value)` 链式调用附加，最终传入函数时通过 `choose_parameter` 提取，若未提供则使用默认值。

这种设计的核心价值是：**向后兼容性**。新增可选参数不会破坏已有调用代码，因为新参数总有默认值。

---

## 三、Property Map：数据与拓扑的解耦

PMP 中几乎所有涉及几何数据的操作都通过 Property Map 进行，而不是直接访问网格的成员变量。

**顶点坐标**通过 `vertex_point_map` 传入：

```cpp
// 默认使用网格内置的点坐标
auto vpm = get(CGAL::vertex_point, mesh);

// 也可以传入外部坐标数组（不修改网格本身）
std::map<vertex_descriptor, Point_3> external_points;
PMP::compute_normals(mesh, CGAL::parameters::vertex_point_map(
    boost::make_assoc_property_map(external_points)));
```

**法向量、曲率等派生属性**通过输出 Property Map 返回：

```cpp
// 将法向量写入 Surface_mesh 的动态属性
auto [vnormals, created] = mesh.add_property_map<vertex_descriptor, Vector_3>("v:normals");
PMP::compute_vertex_normals(mesh, vnormals);
```

这种设计的哲学与 Polygon 包的"无缓存原则"一脉相承：PMP 不在网格内部存储派生属性，而是将存储策略的决定权交给用户。用户可以选择将结果存入网格的动态属性、外部 map，或者直接丢弃（用 `CGAL::Emptyset_iterator`）。

---

## 四、算法的模块化与包的拆分

CGAL 6.2 对 PMP 进行了重大重组，将原本单一的大包拆分为多个专注包：

- `Polygon_mesh_processing`（核心包）：谓词、法向量、曲率、距离、连通分量。
- `PMP_Boolean_operations`：布尔运算、裁剪、切片。
- `PMP_Remeshing`：重网格化、细化、简化、平滑。
- `PMP_Mesh_repair`：网格修复、洞填充、退化消除、边界缝合。

这种拆分体现了**单一职责原则**：每个包只做一件事，依赖关系更清晰，用户只需引入实际需要的功能。从代码层面看，`polygon_mesh_processing.h` 仍然是一个便利头文件，`#include` 所有子模块，保持向后兼容。

---

## 五、Triangle Soup 与 Polygon Mesh 的双轨支持

PMP 同时支持两种输入格式：

- **Polygon Mesh**：有明确拓扑结构的半边网格，面之间共享边和顶点。
- **Triangle Soup**：无拓扑的三角形集合，仅有顶点坐标和面索引数组。

许多算法（如自相交检测 `self_intersections`、法向量计算）对两种格式都有实现。Triangle Soup 的支持通过 `Triangle_mesh_and_triangle_soup_wrapper` 这样的内部适配器实现，对外暴露统一的接口。

这种设计的实用价值在于：从文件读入的原始数据往往是 Triangle Soup 格式，用户可以先用 PMP 的 soup 版本函数做预处理（如 `orient_polygon_soup`、`repair_polygon_soup`），再转换为 Polygon Mesh（`polygon_soup_to_polygon_mesh`）进行后续操作。

---

## 六、精确性与性能的分层策略

PMP 的算法在精确性上采用分层策略：

**谓词层**（如自相交检测、点包含测试）：使用精确谓词，推荐配合 `Exact_predicates_inexact_constructions_kernel`（EPICK）。谓词的正确性由 Kernel 的 filtered 机制保证，性能损失可控。

**构造层**（如重网格化、洞填充）：涉及几何构造（新点坐标计算），通常使用 EPICK 即可满足实际需求。若需要绝对精确的构造结果，可使用 `Exact_predicates_exact_constructions_kernel`（EPECK），但性能代价显著。

**度量层**（如面积、体积、边长计算）：使用浮点运算，结果是近似值。`measure.h` 中的函数明确说明返回值类型为 `FT`（Kernel 的数值类型），精度由所选 Kernel 决定。

这种分层策略与 Polygon 包"精确性由 Kernel 保证"的哲学一致，但 PMP 因为算法复杂度更高，对 Kernel 选择的性能影响更为敏感，因此在文档中对每个函数都明确说明了推荐的 Kernel。

---

## 七、并行化支持

PMP 的部分计算密集型算法支持 TBB（Threading Building Blocks）并行化：

```cpp
// self_intersections.h 中的并行支持
#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#endif
```

并行化通过编译时宏控制，不影响串行版本的接口。这是 Polygon 包所没有的特性，反映了 PMP 处理大规模网格（百万面以上）的实际需求。

---

## 八、与 Polygon 包的设计对比

| 维度 | Polygon 包 | PMP 包 |
|---|---|---|
| 泛型机制 | Traits 参数化 | BGL 图概念 |
| 可选参数 | 固定签名 | Named Parameters |
| 数据存储 | 无缓存，每次重算 | Property Map，用户控制 |
| 维度 | 2D | 3D |
| 并行化 | 无 | TBB 可选 |
| 包规模 | 单一包 | 拆分为多个子包 |

两个包都遵循"算法与数据结构分离"和"精确性由 Kernel 保证"的核心哲学，但 PMP 因为处理对象（三维网格）的复杂性，在 API 设计上引入了更多的灵活性机制。
