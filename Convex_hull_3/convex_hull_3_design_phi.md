# CGAL Convex_hull_3 模板设计与 C++ 思想解析

> 以 Convex_hull_3 为主线，串联 CGAL 中反复出现的模板设计模式。

---

## 1. Traits：Policy-Based Design 的核心

### 1.1 问题的起点

几何算法依赖"判断"（谓词），而判断的精度和实现方式应该与算法逻辑解耦。比如"点是否在平面正侧"这个判断，可以用浮点数快速计算，也可以用精确算术保证正确性。算法本身不应该关心这个细节。

Traits 类就是这个解耦层，它是 C++ **Policy-Based Design** 的直接体现。

### 1.2 Convex_hull_traits_3 的结构

```cpp
// Convex_hull_traits_3.h
template <class R_,
          class Polyhedron = Default,
          class Has_filtered_predicates_tag = Boolean_tag<
              std::is_floating_point<typename R_::FT>::type::value &&
              R_::Has_filtered_predicates_tag::value >>
class Convex_hull_traits_3 {
public:
    typedef R_                    R;
    typedef typename R::Point_3   Point_3;
    typedef Point_triple<R>       Plane_3;   // ← 关键：平面不用 R::Plane_3

    typedef Point_triple_has_on_positive_side_3<R>          Has_on_positive_side_3;
    typedef Point_triple_less_signed_distance_to_plane_3<R> Less_signed_distance_to_plane_3;

    // 工厂方法
    Has_on_positive_side_3 has_on_positive_side_3_object() const { return {}; }
    Less_signed_distance_to_plane_3 less_signed_distance_to_plane_3_object() const { return {}; }
};
```

**三个模板参数各自的职责：**

- `R_`：几何 Kernel，提供点、向量、基础谓词
- `Polyhedron`：输出网格类型（默认 `Polyhedron_3<R>`），用 `Default::Get` 处理默认值
- `Has_filtered_predicates_tag`：**编译期自动检测**是否启用过滤谓词

第三个参数的计算方式值得注意：

```cpp
Boolean_tag<
    std::is_floating_point<typename R_::FT>::type::value &&
    R_::Has_filtered_predicates_tag::value
>
```

这是纯编译期逻辑：如果 Kernel 的数值类型是浮点数，且 Kernel 本身支持过滤谓词，则自动启用三层精度过滤。用户不需要手动指定，Traits 自己推断。

### 1.3 Point_triple：重新定义"平面"

`Convex_hull_traits_3` 把 `Plane_3` 定义为 `Point_triple<R>`（三个点），而不是 Kernel 的 `Plane_3`（ax+by+cz+d=0 形式）。

```cpp
template <class R_>
class Point_triple {
    Point_3 p_, q_, r_;
public:
    const Point_3& p() const { return p_; }
    const Point_3& q() const { return q_; }
    const Point_3& r() const { return r_; }
};
```

原因：QuickHull 算法中，面由三个顶点定义，直接用这三个点表示平面，避免了"构造平面方程"这一步的精度损失。谓词直接用三点做 `Orientation_3` 判断，比先构造平面再判断更精确、更快。

---

## 2. Rebind_TDS：模板的"重绑定"模式

### 2.1 问题

`Triangulation_data_structure_2<Vb, Fb>` 需要 `Vb`（顶点基类）和 `Fb`（面基类）作为模板参数。但 `Vb` 和 `Fb` 内部需要引用 TDS 本身（比如 `Vertex_handle` 的类型依赖 TDS）。这是一个循环依赖。

### 2.2 解决方案：Rebind_TDS

每个 base 类都提供一个内嵌模板 `Rebind_TDS`，允许 TDS 在实例化时"重绑定"自身：

```cpp
// Convex_hull_face_base_2.h
template <typename GT, typename Fb = Triangulation_ds_face_base_2<>>
class Convex_hull_face_base_2 : public Fb {
public:
    // 额外数据：外部点集合 + pending_facets 中的迭代器
    std::list<typename GT::Point_3> points;
    std::list<Face_handle>::iterator it;
    int _info = 0;

    // Rebind_TDS：当 TDS 实例化时，用实际的 TDS2 替换占位符
    template <typename TDS2>
    struct Rebind_TDS {
        typedef typename Fb::template Rebind_TDS<TDS2>::Other Fb2;
        typedef Convex_hull_face_base_2<GT, Fb2>              Other;
    };
};
```

TDS 的实例化过程：

```cpp
// convex_hull_3.h 中
typedef Triangulation_data_structure_2<
    Convex_hull_vertex_base_2<GT3_for_CH3<Traits>>,
    Convex_hull_face_base_2<Traits>
> Tds;

// TDS 内部会做：
typedef typename Vb::template Rebind_TDS<Tds>::Other Vertex_base;
typedef typename Fb::template Rebind_TDS<Tds>::Other Face_base;
```

这个模式解决了"类型需要引用自身所在容器"的循环依赖，是 CGAL 中 TDS 系列（TDS_2、TDS_3）的通用机制。

### 2.3 GT3_for_CH3：类型适配器

```cpp
template <typename GT>
struct GT3_for_CH3 {
    typedef typename GT::Point_3 Point_2;  // 把 Point_3 重命名为 Point_2
};
```

`Convex_hull_vertex_base_2` 期望 Traits 提供 `Point_2`，但 QuickHull 处理的是 3D 点。这个小 struct 做了一个类型重命名，让 3D 的顶点基类能用在 2D 的 TDS 框架里。这是**类型适配器（Type Adapter）**模式的最小化实现。

---

## 3. Forward_functor：谓词的透明转发

### 3.1 场景

`Extreme_points_traits_adapter_3` 处理的是"顶点索引"而非"点坐标"。算法内部调用谓词时传入的是索引，但底层谓词期望坐标。需要一个透明的转发层。

### 3.2 实现

```cpp
// Extreme_points_traits_adapter_3.h
template <class F, class PointPropertyMap>
struct Forward_functor : public F {
    PointPropertyMap vpm_;

    Forward_functor(const PointPropertyMap& vpm, const F& f) : F(f), vpm_(vpm) {}

    template <class Vertex>
    decltype(auto) operator()(const Vertex& p, const Vertex& q) const {
        return static_cast<const F*>(this)->operator()(get(vpm_, p), get(vpm_, q));
    }

    template <class Vertex>
    decltype(auto) operator()(const Vertex& p, const Vertex& q, const Vertex& r) const {
        return static_cast<const F*>(this)->operator()(get(vpm_, p), get(vpm_, q), get(vpm_, r));
    }
    // 4参数版本同理
};
```

**设计要点：**

- 继承自 `F`（原始谓词），保留其所有接口
- 重载 `operator()` 拦截调用，通过 `get(vpm_, p)` 把索引转换为坐标
- `decltype(auto)` 完美转发返回类型，不引入额外拷贝
- `static_cast<const F*>(this)->operator()` 显式调用父类版本，避免无限递归

然后 `Extreme_points_traits_adapter_3` 用这个模板批量生成所有谓词的转发版本：

```cpp
typedef Forward_functor<typename Base_traits::Equal_3,    PointPropertyMap> Equal_3;
typedef Forward_functor<typename Base_traits::Collinear_3, PointPropertyMap> Collinear_3;
typedef Forward_functor<typename Base_traits::Coplanar_3,  PointPropertyMap> Coplanar_3;
// ...
```

这是**装饰器模式（Decorator Pattern）**在 C++ 模板中的实现：不修改原始谓词，在外层包一个转发层。

---

## 4. Is_on_positive_side_of_plane_3：三层精度过滤

这是整个包中最精密的工程实现，集中体现了 CGAL 对"速度与正确性"的权衡。

### 4.1 三个特化版本

```cpp
// 通用版本（非浮点 Kernel）
template <class Traits, class Is_CK = std::false_type>
class Is_on_positive_side_of_plane_3 { ... };

// Extreme_points_traits_adapter_3 的特化（透明转发）
template <class Base_traits, class VPM, class Is_CK>
class Is_on_positive_side_of_plane_3<Extreme_points_traits_adapter_3<VPM,Base_traits>, Is_CK>
    : public Is_on_positive_side_of_plane_3<Base_traits> { ... };

// 浮点 Kernel 的特化（三层过滤）
template <class Kernel, class P>
class Is_on_positive_side_of_plane_3<Convex_hull_traits_3<Kernel, P, Tag_true>, std::true_type> { ... };
```

通过模板特化，编译器自动选择最合适的实现，用户无感知。

### 4.2 三层过滤的实现逻辑

```
构造时（一次性计算）：
    m10 = pqy*prz - pry*pqz   ← 法向量分量（叉积）
    m20 = pqx*prz - prx*pqz
    m21 = pqx*pry - prx*pqy
    Maxx, Maxy, Maxz           ← 坐标绝对值的最大值（用于误差界）

调用时（每次判断）：
    第一层：静态过滤
        det = psx*m10 - m20*psy + m21*psz
        eps = 5.11e-15 * maxx * maxy * maxz   ← 理论误差界
        if |det| > eps → 直接返回（覆盖 ~99% 情况）

    第二层：区间算术（懒惰初始化）
        用 Interval_nt_advanced 重算
        if 区间不含 0 → 确定结果

    第三层：精确算术（懒惰初始化，用指针延迟构造）
        用 Simple_cartesian<Exact_field_selector<double>::Type>
        保证绝对正确
```

**懒惰初始化**的实现：

```cpp
mutable Vector_plus_point<Approx_K> ak_plane;   // 区间算术的平面
mutable Vector_plus_point<Exact_K>* ek_plane_ptr; // 精确算术的平面（指针，按需构造）

// 用 infinity() 作为"未初始化"的哨兵值
ak_plane.vector = Vector_3(Interval_nt_advanced(0., std::numeric_limits<double>::infinity()), 0., 0.);

// 调用时检查
if (ak_plane.vector.x().sup() == std::numeric_limits<double>::infinity()) {
    // 第一次需要区间算术时才构造
    ak_plane.vector = cross_product(to_AK(q)-to_AK(p), to_AK(r)-to_AK(p));
}
```

`mutable` 关键字允许在 `const` 的 `operator()` 中修改缓存状态，这是 C++ 中实现**逻辑常量性（logical constness）**的标准手法。

---

## 5. 模板特化控制算法路径

### 5.1 BOOST_MPL_HAS_XXX_TRAIT

```cpp
// convex_hull_3.h
BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_type_Collinear_3, Collinear_3, false)
```

这个宏生成一个 trait 类，检测某个类型是否有名为 `Collinear_3` 的嵌套类型。用于在编译期判断 Traits 是否提供了某个可选谓词，从而选择不同的算法路径。

### 5.2 if constexpr 在新代码中的使用

`helpers.h` 中的 `extreme_point_3_wrapper` 展示了现代 C++ 风格：

```cpp
template <class Convex, class Direction_3, class NamedParameters>
auto extreme_point_3_wrapper(const Convex& c, const Direction_3& dir, const NamedParameters& np) {
    if constexpr (is_instance_of_CHH<Convex>) {
        // Convex_hull_hierarchy_3 的特殊路径
        return c.extreme_point_3(dir);
    } else if constexpr (CGAL::IO::internal::is_Range_v<Convex>) {
        // 点集范围的路径
        // ...
    } else {
        // 普通 FaceGraph 的路径
        // ...
    }
}
```

同一个函数处理三种完全不同的输入类型，编译期分支，零运行时开销。与旧代码中的 tag dispatch 相比，`if constexpr` 更直观，但可扩展性稍弱（第三方无法注入新分支）。

---

## 6. Named Parameters：可选参数的现代处理

CGAL 5.x 引入了 Named Parameters 机制，替代了大量重载函数：

```cpp
// 用户调用
CGAL::convex_hull_3(points.begin(), points.end(), mesh,
    CGAL::parameters::geom_traits(my_traits));

// 内部提取
using CGAL::parameters::choose_parameter;
using CGAL::parameters::get_parameter;

auto traits = choose_parameter(get_parameter(np, internal_np::geom_traits),
                               Default_traits{});
```

**实现原理**：Named Parameters 是一个链式的类型列表，每个参数是一个 `(tag, value)` 对。`get_parameter` 在编译期遍历这个列表查找对应 tag，`choose_parameter` 在找不到时返回默认值。整个过程在编译期完成，运行时零开销。

这解决了 C++ 没有默认命名参数的问题，同时保持了类型安全。

---

## 7. 整个包的模板层次总览

```
用户调用
convex_hull_3(first, beyond, P, np)
        │
        ▼
Named Parameters 解析
→ 提取 geom_traits（或使用默认 Convex_hull_traits_3<Kernel>）
        │
        ▼
Convex_hull_traits_3<R, Polyhedron, Has_filtered_predicates_tag>
    ├── Point_3 = R::Point_3
    ├── Plane_3 = Point_triple<R>          ← 重定义平面类型
    ├── Has_on_positive_side_3             ← 特化选择精度策略
    └── Less_signed_distance_to_plane_3
        │
        ▼
ch_quickhull_face_graph<InputIterator, PolygonMesh, Traits>
    │
    ├── Triangulation_data_structure_2<
    │       Convex_hull_vertex_base_2<GT3_for_CH3<Traits>>,  ← 类型适配
    │       Convex_hull_face_base_2<Traits>                  ← 注入算法数据
    │   >
    │       └── Rebind_TDS 解决循环依赖
    │
    ├── Is_on_positive_side_of_plane_3<Traits>
    │       └── 模板特化 → 三层精度过滤 / 通用版本 / 转发版本
    │
    └── copy_face_graph(tds, P)            ← BGL concept 统一输出
            └── 适配 Surface_mesh / Polyhedron_3 / Indexed_triangle_set
```

---

## 8. C++ 设计思想总结

| 模式 | 在包中的体现 | C++ 机制 |
|------|------------|---------|
| Policy-Based Design | Traits 类解耦算法与谓词 | 模板参数 |
| 编译期自动选择策略 | `Has_filtered_predicates_tag` 自动推断 | `Boolean_tag` + 模板特化 |
| 类型重绑定 | `Rebind_TDS` 解决循环依赖 | 内嵌模板 |
| 类型适配器 | `GT3_for_CH3` 重命名 Point_3→Point_2 | typedef |
| 装饰器模式 | `Forward_functor` 透明转发谓词 | 继承 + `operator()` 重载 |
| 逻辑常量性 | 谓词缓存区间/精确算术结果 | `mutable` |
| 懒惰初始化 | 精确算术平面按需构造 | 指针 + 哨兵值 |
| 零开销抽象 | 后置条件编译期可关闭 | 预处理宏 |
| 编译期分支 | `extreme_point_3_wrapper` | `if constexpr` |
| 可扩展分发 | 迭代器类别选择算法 | Tag dispatch |
| 命名参数 | `CGAL::parameters::geom_traits(...)` | 链式类型列表 |
