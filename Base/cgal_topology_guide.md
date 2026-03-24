# CGAL 拓扑数据结构指南

> 理解 CGAL 中点、线、面的存储与遍历方式，是读懂大多数几何包的前提。

---

## 目录

1. [核心概念：三种拓扑表示](#1-核心概念三种拓扑表示)
2. [半边数据结构（HalfedgeDS）](#2-半边数据结构halfedgeds)
3. [Polyhedron_3：经典多面体](#3-polyhedron_3经典多面体)
4. [Surface_mesh：现代网格](#4-surface_mesh现代网格)
5. [Triangulation_data_structure_2：三角剖分专用结构](#5-triangulation_data_structure_2三角剖分专用结构)
6. [Circulator：CGAL 的环形迭代器](#6-circulator-cgal-的环形迭代器)
7. [BGL 接口层：统一访问所有网格](#7-bgl-接口层统一访问所有网格)
8. [Property Map：给拓扑元素附加数据](#8-property-map给拓扑元素附加数据)
9. [选哪个数据结构](#9-选哪个数据结构)

---

## 1. 核心概念：三种拓扑表示

CGAL 中的网格/多面体有三种主要的底层表示，理解它们的区别是一切的起点：

```
HalfedgeDS（半边数据结构）
  └── Polyhedron_3（基于 HalfedgeDS 的多面体，老式 API）

Surface_mesh（基于索引的现代网格，推荐使用）

Triangulation_data_structure_2/3（三角剖分专用，用于 Delaunay/Voronoi）
```

它们都通过 **BGL（Boost Graph Library）concept** 统一暴露给上层算法，所以 `Polygon_mesh_processing` 等包的函数可以同时接受这三种类型。

---

## 2. 半边数据结构（HalfedgeDS）

### 2.1 半边的核心思想

半边结构是 CGAL 中最基础的拓扑表示。每条几何边被拆分为**两条方向相反的半边**：

```
顶点 A ──────────────→ 顶点 B
         半边 h
顶点 A ←────────────── 顶点 B
         半边 h->opposite()
```

每条半边存储四个指针：

```
halfedge h:
  ├── opposite()  → 对面的半边（同一条边，反方向）
  ├── next()      → 同一个面内的下一条半边（逆时针）
  ├── prev()      → 同一个面内的上一条半边（可选）
  ├── vertex()    → 半边指向的顶点（终点）
  └── face()      → 半边所属的面（可选）
```

顶点存储：
```
vertex v:
  └── halfedge() → 以该顶点为起点的任意一条半边
```

面存储：
```
face f:
  └── halfedge() → 该面边界上的任意一条半边
```

### 2.2 HalfedgeDS_halfedge_base 的模板参数

```cpp
// 来自 HalfedgeDS_halfedge_base.h
template < class Refs,
           class TP = Tag_true,   // 是否支持 prev()
           class TV = Tag_true,   // 是否支持 vertex()
           class TF = Tag_true>   // 是否支持 face()
class HalfedgeDS_halfedge_base;
```

这三个 `Tag_true/false` 控制半边存储哪些指针。最小化配置（全 false）只存 `opposite` 和 `next`，节省内存。完整配置（全 true）支持所有遍历操作。

### 2.3 HalfedgeDS_vertex_base 的模板参数

```cpp
// 来自 HalfedgeDS_vertex_base.h
template < class Refs,
           class T = Tag_true,    // 是否支持 halfedge()（顶点→半边的链接）
           class P = Tag_false>   // 点坐标类型，Tag_false 表示不存坐标
class HalfedgeDS_vertex_base;
```

典型用法：`HalfedgeDS_vertex_base<Refs, Tag_true, Point_3>` 表示顶点存一条半边引用和一个 3D 坐标。

### 2.4 遍历模式

**绕顶点一圈的所有半边**（star of a vertex）：

```cpp
// 从顶点 v 出发，绕一圈
Halfedge_handle h = v->halfedge();
do {
    // 处理 h
    h = h->next()->opposite();  // 跳到下一条以 v 为起点的半边
} while (h != v->halfedge());
```

**绕面一圈的所有半边**：

```cpp
Halfedge_handle h = f->halfedge();
do {
    // 处理 h->vertex()（面的顶点）
    h = h->next();
} while (h != f->halfedge());
```

---

## 3. Polyhedron_3：经典多面体

### 3.1 定义

```cpp
#include <CGAL/Polyhedron_3.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron;
```

`Polyhedron_3` 是对 `HalfedgeDS` 的高层封装，提供了更友好的 API。

### 3.2 关键类型

```cpp
Polyhedron::Vertex_handle        // 顶点句柄（指针语义）
Polyhedron::Halfedge_handle      // 半边句柄
Polyhedron::Facet_handle         // 面句柄
Polyhedron::Vertex_iterator      // 顶点迭代器
Polyhedron::Halfedge_iterator    // 半边迭代器
Polyhedron::Facet_iterator       // 面迭代器
Polyhedron::Halfedge_around_vertex_circulator   // 绕顶点的环形迭代器
Polyhedron::Halfedge_around_facet_circulator    // 绕面的环形迭代器
```

### 3.3 遍历示例

```cpp
Polyhedron P;
// ... 构建 P

// 遍历所有顶点
for (auto v = P.vertices_begin(); v != P.vertices_end(); ++v) {
    K::Point_3 pt = v->point();
}

// 遍历面的所有顶点（用 circulator）
for (auto f = P.facets_begin(); f != P.facets_end(); ++f) {
    auto hc = f->facet_begin();  // Halfedge_around_facet_circulator
    do {
        K::Point_3 pt = hc->vertex()->point();
        ++hc;
    } while (hc != f->facet_begin());
}

// 遍历顶点的所有邻居
Polyhedron::Vertex_handle v = ...;
auto vc = v->vertex_begin();  // Halfedge_around_vertex_circulator
do {
    Polyhedron::Vertex_handle neighbor = vc->opposite()->vertex();
    ++vc;
} while (vc != v->vertex_begin());
```

### 3.4 Polyhedron_3 的局限

- 句柄是指针，删除元素后句柄失效
- 不支持直接附加自定义属性（需要继承 Items 类）
- 不支持非流形拓扑（每条边恰好属于两个面）
- 现代代码推荐用 `Surface_mesh` 替代

---

## 4. Surface_mesh：现代网格

### 4.1 设计哲学：索引而非指针

`Surface_mesh` 的核心设计与 `Polyhedron_3` 完全不同：用**整数索引**代替指针句柄。

```cpp
// 来自 Surface_mesh.h
class SM_Index<T> {
    uint32_t idx_;  // 底层就是一个整数
};

class SM_Vertex_index   : public SM_Index<SM_Vertex_index>   {};
class SM_Halfedge_index : public SM_Index<SM_Halfedge_index> {};
class SM_Edge_index     : public SM_Index<SM_Edge_index>     {};
class SM_Face_index     : public SM_Index<SM_Face_index>     {};
```

索引的优势：
- 可以存入数组、哈希表，序列化
- 删除元素后索引仍然有效（元素被标记为 deleted，索引不变）
- 支持 `garbage_collection()` 压缩存储

### 4.2 内部存储结构

所有拓扑数据存在连续的 `std::vector` 中：

```
vconn_[i]  → 顶点 i 的连接信息（出发半边的索引）
hconn_[i]  → 半边 i 的连接信息（face, next, prev, vertex）
fconn_[i]  → 面 i 的连接信息（起始半边的索引）
```

这种 SoA（Structure of Arrays）布局对缓存友好，遍历性能优于 `Polyhedron_3` 的指针链表。

### 4.3 关键类型

```cpp
typedef CGAL::Surface_mesh<Point_3> Mesh;

Mesh::Vertex_index    v;   // 顶点索引
Mesh::Halfedge_index  h;   // 半边索引
Mesh::Edge_index      e;   // 边索引（= 一对半边）
Mesh::Face_index      f;   // 面索引
```

注意：`Surface_mesh` 显式区分了 `Edge_index` 和 `Halfedge_index`，而 `Polyhedron_3` 没有独立的 Edge 类型。

### 4.4 遍历示例

```cpp
Mesh sm;
// ... 构建 sm

// 遍历所有顶点
for (Mesh::Vertex_index v : sm.vertices()) {
    Point_3 pt = sm.point(v);
}

// 遍历面的所有顶点
for (Mesh::Face_index f : sm.faces()) {
    // halfedges_around_face 返回一个 range
    for (Mesh::Halfedge_index h : sm.halfedges_around_face(sm.halfedge(f))) {
        Mesh::Vertex_index v = sm.target(h);
        Point_3 pt = sm.point(v);
    }
}

// 遍历顶点的所有邻居
Mesh::Vertex_index v = ...;
for (Mesh::Vertex_index nb : sm.vertices_around_target(sm.halfedge(v))) {
    // nb 是 v 的邻居顶点
}
```

### 4.5 半边操作 API

```cpp
sm.halfedge(v)          // 顶点 v 的一条出发半边
sm.halfedge(f)          // 面 f 的一条边界半边
sm.opposite(h)          // h 的对边
sm.next(h)              // 同面内的下一条半边
sm.prev(h)              // 同面内的上一条半边
sm.target(h)            // h 的终点顶点
sm.source(h)            // h 的起点顶点（= target(opposite(h))）
sm.face(h)              // h 所属的面（border 半边返回 null_face()）
sm.is_border(h)         // h 是否是边界半边（无面）
sm.edge(h)              // h 对应的无向边
sm.halfedge(e, 0/1)     // 边 e 的两条半边
```

---

## 5. Triangulation_data_structure_2：三角剖分专用结构

### 5.1 与网格结构的区别

`TDS_2` 不是通用网格，专门为三角剖分设计：
- 所有面都是三角形（3 个顶点，3 条半边）
- 有一个特殊的**无限顶点**（infinite vertex），用于处理凸包边界
- 面通过**顶点索引**（0/1/2）而非半边来访问邻居

```cpp
// 来自 Triangulation_data_structure_2.h
template < class Vb = Triangulation_ds_vertex_base_2<>,
           class Fb = Triangulation_ds_face_base_2<> >
class Triangulation_data_structure_2;
```

### 5.2 面的邻居访问

TDS_2 的面存储三个邻居面的指针，通过**对面顶点的索引**来标识：

```
面 f 有顶点 v0, v1, v2（逆时针）
f->neighbor(i) = 与 f 共享对面顶点 vi 的对边的邻居面
f->vertex(i)   = 第 i 个顶点
f->index(v)    = 顶点 v 在面 f 中的索引（0/1/2）
```

这与半边结构完全不同——TDS_2 没有显式的半边对象，边由 `(face, index)` 对表示。

### 5.3 遍历示例

```cpp
typedef CGAL::Triangulation_2<K> Triangulation;
Triangulation T;
// ... 插入点

// 遍历所有有限面
for (auto f = T.finite_faces_begin(); f != T.finite_faces_end(); ++f) {
    for (int i = 0; i < 3; ++i) {
        K::Point_2 pt = f->vertex(i)->point();
    }
}

// 遍历顶点的所有邻居面（circulator）
Triangulation::Vertex_handle v = ...;
auto fc = T.incident_faces(v);  // Face_circulator
auto done = fc;
do {
    if (!T.is_infinite(fc)) {
        // 处理有限面
    }
    ++fc;
} while (fc != done);

// 遍历顶点的所有邻居顶点
auto vc = T.incident_vertices(v);  // Vertex_circulator
auto vdone = vc;
do {
    if (!T.is_infinite(vc)) {
        K::Point_2 pt = vc->point();
    }
    ++vc;
} while (vc != vdone);
```

### 5.4 无限顶点与无限面

TDS_2 在凸包外围添加一个虚拟的无限顶点，使得每条边界边都有两个邻居面（一个有限，一个无限）。这简化了边界处理：

```cpp
T.infinite_vertex()          // 获取无限顶点
T.is_infinite(v)             // 判断顶点是否是无限顶点
T.is_infinite(f)             // 判断面是否包含无限顶点
T.finite_faces_begin/end()   // 只遍历有限面
T.finite_vertices_begin/end() // 只遍历有限顶点
```

---

## 6. Circulator：CGAL 的环形迭代器

### 6.1 为什么需要 Circulator

标准迭代器有 `begin()` 和 `end()`，适合线性序列。但拓扑遍历（绕顶点一圈、绕面一圈）是**循环的**，没有自然的起点和终点。CGAL 的 Circulator 解决这个问题。

### 6.2 Circulator 的使用模式

```cpp
// 标准 do-while 模式（不能用 for 循环，因为没有 end）
auto c = get_some_circulator();
if (c != nullptr) {  // 空检查
    do {
        // 处理 *c
        ++c;
    } while (c != get_some_circulator());  // 回到起点时停止
}
```

注意：`c != nullptr` 检查是必要的，空的 circulator（如孤立顶点）会导致无限循环。

### 6.3 Circulator 转 Iterator Range

有时需要把 circulator 用在 range-based for 或 STL 算法中：

```cpp
#include <CGAL/circulator.h>

auto c = f->facet_begin();
// 转换为 [begin, end) 的 iterator range
auto range = CGAL::make_range(c);  // 或用 Container_from_circulator
for (auto& h : range) { ... }
```

---

## 7. BGL 接口层：统一访问所有网格

### 7.1 核心思想

CGAL 通过特化 `boost::graph_traits` 让所有网格类型满足 BGL 的 graph concept，从而让算法与具体数据结构解耦。

```cpp
// 对 Surface_mesh 的特化（graph_traits_Surface_mesh.h）
namespace boost {
template <class P>
struct graph_traits<CGAL::Surface_mesh<P>> {
    typedef SM_Vertex_index   vertex_descriptor;
    typedef SM_Halfedge_index halfedge_descriptor;
    typedef SM_Edge_index     edge_descriptor;
    typedef SM_Face_index     face_descriptor;
    // ...
};
}
```

同样的特化存在于 `Polyhedron_3`、`HalfedgeDS`、`Linear_cell_complex` 等。

### 7.2 BGL 的自由函数 API

通过 BGL 接口，所有网格类型使用相同的自由函数：

```cpp
// 适用于任何满足 FaceGraph concept 的类型 G
template <class G>
void example(const G& g) {
    // 遍历顶点
    for (auto v : vertices(g)) { ... }

    // 遍历面
    for (auto f : faces(g)) { ... }

    // 半边操作
    auto h = halfedge(f, g);
    auto h_opp = opposite(h, g);
    auto h_next = next(h, g);
    auto v = target(h, g);

    // 获取点坐标（通过 property map）
    auto vpm = get(CGAL::vertex_point, g);
    auto pt = get(vpm, v);
}
```

### 7.3 常用 BGL Concept 层次

```
Graph
  └── HalfedgeGraph        → 有半边操作（opposite, next, prev）
        └── FaceGraph       → 有面操作（face, halfedge(f,g)）
              └── MutableFaceGraph  → 可修改（add_vertex, add_face 等）
```

读 CGAL 算法的模板参数时，看到 `class TriangleMesh` 或 `class PolygonMesh`，通常要求满足 `FaceGraph` 或 `MutableFaceGraph`。

### 7.4 BGL 迭代器范围

```cpp
// CGAL/boost/graph/iterator.h 提供了大量便利的 range
vertices(g)                          // 所有顶点
halfedges(g)                         // 所有半边
edges(g)                             // 所有边
faces(g)                             // 所有面
halfedges_around_face(h, g)          // 绕面 h 所在面的所有半边
halfedges_around_target(h, g)        // 以 target(h) 为终点的所有半边
vertices_around_target(h, g)         // target(h) 的所有邻居顶点
faces_around_target(h, g)            // target(h) 的所有邻居面
```

---

## 8. Property Map：给拓扑元素附加数据

### 8.1 概念

Property Map 是 CGAL 中给顶点/半边/面附加任意数据的标准机制，解耦了"拓扑结构"和"附加属性"。

### 8.2 Surface_mesh 的内置 Property Map

```cpp
Mesh sm;

// 内置的点坐标 property map（每个网格都有）
auto vpm = sm.points();  // 或 get(CGAL::vertex_point, sm)
Point_3 pt = get(vpm, v);
put(vpm, v, new_point);

// 添加自定义属性
auto [vnormals, created] = sm.add_property_map<Mesh::Vertex_index, Vector_3>("v:normal");
// 使用
Vector_3 n = get(vnormals, v);
put(vnormals, v, new_normal);

// 删除属性
sm.remove_property_map(vnormals);
```

### 8.3 Property Map 的内部实现

来自 `Surface_mesh/Properties.h`：

```cpp
// 所有属性数组的基类（类型擦除）
class Base_property_array {
    virtual void reserve(size_t n) = 0;
    virtual void resize(size_t n) = 0;
    virtual void push_back() = 0;
    // ...
};

// 具体类型的属性数组
template <class T>
class Property_array : public Base_property_array {
    std::vector<T> data_;  // 连续存储，索引对应网格元素
};
```

属性数组与网格元素数组等长，通过相同的整数索引访问。当网格增删元素时，所有属性数组同步 resize。

### 8.4 外部 Property Map（不修改网格）

有时不想修改网格结构，可以用外部 map：

```cpp
// 用 std::unordered_map 作为 property map
std::unordered_map<Mesh::Vertex_index, double> dist_map;
auto dist_pmap = boost::make_assoc_property_map(dist_map);

// 用 std::vector 作为 property map（索引必须连续）
std::vector<int> labels(sm.num_vertices());
auto label_pmap = CGAL::make_property_map(labels);
```

---

## 9. 选哪个数据结构

| 场景 | 推荐 | 原因 |
|------|------|------|
| 通用网格处理（PMP 等算法） | `Surface_mesh` | 现代 API，性能好，property map 支持完善 |
| 需要自定义顶点/面数据类型 | `Polyhedron_3` + 自定义 Items | 通过继承 Items 类扩展 |
| Delaunay 三角化 / Voronoi | `Triangulation_2/3` | 专用结构，支持无限顶点，插入删除高效 |
| 只读算法，输入来自外部 | 任意满足 `FaceGraph` 的类型 | BGL 接口统一 |
| 需要动态插入删除点 | `Delaunay_triangulation_3` | 基于 TDS_3，支持增量更新 |
| 非流形拓扑（如点云） | `Linear_cell_complex` | 支持任意维度和非流形 |

### 关键判断问题

1. **需要修改拓扑吗？** → 需要 `MutableFaceGraph`，`Surface_mesh` 或 `Polyhedron_3`
2. **需要附加自定义属性吗？** → `Surface_mesh` 的 property map 最方便
3. **输入是三角剖分/Voronoi 吗？** → 用 `Triangulation_*` 系列
4. **只是调用 PMP/AABB 等算法？** → 用 `Surface_mesh`，最通用

---

## 附录：常见遍历模式速查

### Surface_mesh

```cpp
// 面的顶点
for (auto h : sm.halfedges_around_face(sm.halfedge(f)))
    auto v = sm.target(h);

// 顶点的邻居顶点
for (auto v2 : sm.vertices_around_target(sm.halfedge(v)))
    ...;

// 顶点的邻居面
for (auto f2 : sm.faces_around_target(sm.halfedge(v)))
    ...;
```

### Polyhedron_3

```cpp
// 面的顶点（circulator）
auto hc = f->facet_begin();
do { auto pt = hc->vertex()->point(); } while (++hc != f->facet_begin());

// 顶点的邻居（circulator）
auto vc = v->vertex_begin();
do { auto nb = vc->opposite()->vertex(); } while (++vc != v->vertex_begin());
```

### Triangulation_2

```cpp
// 顶点的邻居面（circulator）
auto fc = T.incident_faces(v), done = fc;
do { if (!T.is_infinite(fc)) { ... } } while (++fc != done);

// 顶点的邻居顶点（circulator）
auto vc = T.incident_vertices(v), done = vc;
do { if (!T.is_infinite(vc)) { ... } } while (++vc != done);
```
