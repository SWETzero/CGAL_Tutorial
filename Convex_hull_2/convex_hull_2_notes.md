# CGAL Convex_hull_2 包深度解析

> 基于 CGAL main 分支源码分析，聚焦算法设计、C++ 现代模板实践与工程取舍。

---

## 目录

1. [包的整体结构](#1-包的整体结构)
2. [Traits 体系：泛型几何的基础](#2-traits-体系泛型几何的基础)
3. [迭代器分发：基于 tag dispatch 的算法选择](#3-迭代器分发基于-tag-dispatch-的算法选择)
4. [算法详解](#4-算法详解)
5. [算法选取的工程逻辑](#5-算法选取的工程逻辑)
6. [后置条件与正确性验证](#6-后置条件与正确性验证)
7. [Benchmark 与性能测试](#7-benchmark-与性能测试)
8. [现代 C++ 模板实践总结](#8-现代-c-模板实践总结)

---

## 1. 包的整体结构

```
Convex_hull_2/
├── include/CGAL/
│   ├── convex_hull_2.h              ← 主入口，算法分发逻辑
│   ├── ch_akl_toussaint.h           ← 声明
│   ├── ch_bykat.h
│   ├── ch_eddy.h
│   ├── ch_jarvis.h
│   ├── ch_melkman.h
│   ├── ch_graham_andrew.h
│   ├── convex_hull_traits_2.h       ← Traits 定义
│   ├── convexity_check_2.h          ← 后置条件验证
│   └── Convex_hull_2/
│       ├── ch_akl_toussaint_impl.h  ← 实现（与声明分离）
│       ├── ch_bykat_impl.h
│       ├── ch_eddy_impl.h
│       ├── ch_jarvis_impl.h
│       ├── ch_melkman_impl.h
│       └── ch_graham_andrew_impl.h
└── examples/Convex_hull_2/
    ├── ch_timing.cpp                ← 所有算法的 benchmark
    └── ...
```

**声明与实现分离**的原因：头文件只暴露函数签名，`_impl.h` 包含模板实现体。用户的 `.h` 文件 `#include` 声明头，声明头末尾再 `#include` 对应的 `_impl.h`，这样既保持了接口清晰，又满足了模板必须在头文件中定义的要求。

---

## 2. Traits 体系：泛型几何的基础

### 2.1 设计动机

CGAL 的核心问题：几何算法依赖**谓词**（如"三点是否左转"），而谓词的精度和实现方式应该与算法逻辑解耦。Traits 类就是这个解耦层。

### 2.2 Convex_hull_traits_2 的实现

```cpp
// convex_hull_traits_2.h
template <class K_>
class Convex_hull_traits_2 : public K_
{
public:
  Convex_hull_traits_2() { }
  Convex_hull_traits_2(const K_& k) : K_(k) { }
};
```

这里用了一个极简的设计：**直接继承 Kernel**。CGAL 的 Kernel（如 `Exact_predicates_inexact_constructions_kernel`）本身已经提供了所有需要的谓词对象工厂方法，所以 Traits 不需要额外定义任何东西，继承即可。

### 2.3 各算法依赖的谓词

不同算法需要不同的谓词，这决定了它们对 Traits 的要求：

| 算法 | 必需谓词 | 说明 |
|------|---------|------|
| Graham-Andrew | `Left_turn_2`, `Less_xy_2`, `Equal_2` | 基础谓词 |
| Akl-Toussaint | 以上 + `Less_yx_2` | 需要找 N/S 极值点 |
| Bykat / Eddy | 以上 + `Compare_signed_distance_to_line_2` | 需要找最远点 |
| Jarvis | `Less_rotate_ccw_2`, `Equal_2`, `Less_xy_2` | 完全不同的谓词集合 |
| Melkman | `Left_turn_2`, `Equal_2` | 最简单 |

`Less_rotate_ccw_2` 是 Jarvis 独有的：给定基准点 `p`，比较 `q1` 和 `q2` 哪个在 `p` 的逆时针方向更靠前。这个谓词在其他算法里没有用到，这也是 Jarvis 无法被其他算法替代的原因之一。

### 2.4 谓词对象工厂模式

CGAL 不直接调用全局函数，而是通过 Traits 获取**谓词函数对象**：

```cpp
// 典型用法（来自 ch_bykat_impl.h）
Left_turn_2 left_turn = ch_traits.left_turn_2_object();
Less_xy_2   less_xy   = ch_traits.less_xy_2_object();

// 然后像函数一样调用
if (left_turn(a, b, c)) { ... }
```

这个模式的好处：谓词对象可以携带状态（比如精度控制参数），而全局函数做不到。

---

## 3. 迭代器分发：基于 tag dispatch 的算法选择

### 3.1 核心分发逻辑

`convex_hull_2.h` 中的分发是整个包最重要的设计决策：

```cpp
// 针对不同迭代器类别，调用不同算法
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
CGAL_convex_hull_points_2(InputIterator first, InputIterator last,
                          OutputIterator result, const Traits& ch_traits,
                          std::input_iterator_tag)          // ← tag 参数
{ return ch_bykat(first, last, result, ch_traits); }

template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
CGAL_convex_hull_points_2(InputIterator first, InputIterator last,
                          OutputIterator result, const Traits& ch_traits,
                          std::forward_iterator_tag)        // ← 不同 tag
{ return ch_akl_toussaint(first, last, result, ch_traits); }

// bidirectional 和 random_access 同样调用 akl_toussaint
```

公开接口提取 tag 并转发：

```cpp
template <class InputIterator, class OutputIterator, class Traits>
OutputIterator
convex_hull_points_2(InputIterator first, InputIterator last,
                     OutputIterator result, const Traits& ch_traits)
{
    typedef std::iterator_traits<InputIterator>::iterator_category Category;
    return CGAL_convex_hull_points_2(first, last, result, ch_traits,
                                     Category());  // ← 构造 tag 对象
}
```

### 3.2 为什么用 tag dispatch 而不是 if constexpr

这段代码写于 C++98/03 时代，tag dispatch 是当时实现编译期分支的标准手段。现代 C++ 可以用 `if constexpr` + `std::is_same_v` 改写，但 tag dispatch 有一个优势：**可扩展性**。第三方可以为自定义迭代器类别特化分发，而 `if constexpr` 链是封闭的。

### 3.3 迭代器类别与算法能力的对应关系

```
input_iterator
  └─ 只能单次遍历，不能回头
  └─ → ch_bykat（内部 std::copy 到 vector，然后多次遍历）

forward_iterator
  └─ 可多次遍历，但只能向前
  └─ → ch_akl_toussaint（需要两次遍历：找极值点 + 分配区域）

bidirectional / random_access
  └─ 同上，akl_toussaint 对 random_access 有额外优化
  └─ → ch_akl_toussaint
```

---

## 4. 算法详解

### 4.1 Graham-Andrew（Andrew's Monotone Chain）

**角色**：基础构件，被 akl_toussaint 调用，也可独立使用。

**核心思路**：

1. 按 x 坐标排序（x 相同按 y）
2. 从左到右扫描构建下凸壳（维护一个栈，遇到右转则弹栈）
3. 从右到左扫描构建上凸壳
4. 合并

**复杂度**：O(n log n)，排序主导。

**实现细节**（来自 `ch_graham_andrew_impl.h`）：

```cpp
// 用迭代器的 vector 作为栈，避免拷贝点数据
std::vector<BidirectionalIterator> S;
// 核心循环：维护左转不变量
while (!left_turn(*beta, *alpha, *iter)) {
    S.pop_back();
    alpha = beta;
    beta = *++stack_rev_iter;
}
```

注意栈里存的是**迭代器**而不是点本身，节省拷贝开销。

**要求**：`BidirectionalIterator`（需要从两端扫描）。

细节：  
 **共线点处理 (Collinear Points)** 在 ``ch_graham_andrew_impl.h`` 中，判定条件是 ``!left_turn`` 还是 ``right_turn``这决定了凸包边上是否保留共线点。``is_ccw_strongly_convex_2`` 验证的是强凸性（即不允许三点共线）。  
 **重复点处理 (Duplicate Points)** 算法如何通过 ``Equal_2`` 谓词跳过重叠点，防止在 ``ch_jarvis`` 旋转判定中出现除零或无限循环。

---

### 4.2 Akl-Toussaint

**角色**：`convex_hull_2()` 对 forward/bidirectional/random_access 迭代器的默认算法。

**核心思路**：预处理过滤 + Graham-Andrew。

**第一步：找四个极值点（O(n)）**

```
N（最北，y 最大）
W（最西，x 最小）  ←——→  E（最东，x 最大）
S（最南，y 最小）
```

这四点构成一个八边形（实际上是四边形，但可以扩展到八边形）。

**第二步：过滤内部点（O(n)）**

将点集分配到四个区域（NE、SE、SW、NW 象限），八边形内部的点直接丢弃。对于均匀随机分布，这一步能过滤掉约 80-90% 的点。

**第三步：对每个区域跑 Graham-Andrew**

```cpp
// 来自 ch_akl_toussaint_impl.h 的结构
ch_nswe_point_with_order(first, last, n, s, w, e, ch_traits);  // 找极值点
ch_akl_toussaint_assign_points_to_regions(...);                 // 分配区域
// 对 region1..4 分别调用 ch_graham_andrew_scan
```

**对 random_access 迭代器的特化优化**：

```cpp
// forward_iterator 版本：需要遍历一次记录位置
template <...> ch_nswe_point_with_order(..., std::forward_iterator_tag)
// 用 pair<position, index> 追踪极值点位置，最后排序

// random_access 版本：直接比较迭代器大小
template <...> ch_nswe_point_with_order(..., std::random_access_iterator_tag)
// 直接对迭代器数组排序，因为 random_access 迭代器支持 <
```

这是 CGAL 中 **iterator tag 特化**的典型用法：同一逻辑，针对不同迭代器能力有不同的高效实现。

**复杂度**：O(n log h)，h 为凸包顶点数；最坏 O(n log n)。

---

### 4.3 Bykat

**角色**：`convex_hull_2()` 对 input_iterator 的默认算法。

**核心思路**：QuickHull 变体，递归分治。

**算法流程**：

1. 找最左点 `a` 和最右点 `b`（需要先把数据复制到 vector）
2. 用直线 `ab` 将点集分为上下两半
3. 在每半中找距离 `ab` 最远的点 `c`
4. 三角形 `abc` 内部的点不可能是凸包顶点，丢弃
5. 递归处理 `ac` 和 `cb` 两段

**实现中的现代 C++ 用法**（来自 `ch_bykat_impl.h`）：

```cpp
// lambda 作为 std::partition 的谓词
l = std::partition(P.begin(), P.end(),
    [&left_turn, &a, &b](const Point_2& p) {
        return left_turn(a, b, p);
    });

// lambda 作为 std::min_element 的比较器（找最远点）
auto less_dist = [&a, &b, &cmp_dist, &less_xy](const Point_2& p1, const Point_2& p2) -> bool {
    CGAL::Comparison_result res = cmp_dist(a, b, p1, p2);
    if (res == CGAL::EQUAL) return less_xy(p1, p2);
    return (res == CGAL::SMALLER);
};
Point_2 c = *std::min_element(l, r, less_dist);
```

用 `std::partition` 原地分割，避免额外内存分配；用 lambda 捕获 traits 对象，保持泛型性。

**`ch_bykat_with_threshold` 变体**：当子问题规模小于阈值（`CGAL_ch_THRESHOLD = 10`）时，切换到 Graham-Andrew，避免递归开销。

**复杂度**：平均 O(n log h)，最坏 O(n²)（所有点共线或退化输入）。

---

### 4.4 Eddy

**角色**：独立提供，不作为默认算法。

**与 Bykat 的关系**：思路相同（QuickHull 递归分治），Eddy 是更早的实现。主要区别：

- Eddy 用 `std::list` 存储中间结果，Bykat 用 `std::vector` + `std::partition` 原地操作
- Bykat 的内存访问模式更友好（vector 连续内存 vs list 链表跳转）
- Bykat 有 `with_threshold` 变体，Eddy 没有

Eddy 保留在库中的原因：历史遗留 + 作为参考实现，方便对比测试。

**复杂度**：与 Bykat 相同，平均 O(n log h)，最坏 O(n²)。

---

### 4.5 Jarvis March（Gift Wrapping）

**角色**：独立提供，不作为默认算法，但有独特的接口价值。

**核心思路**：从最左点出发，每次用 `std::min_element` 找角度最小的下一个点。

**实现**（来自 `ch_jarvis_impl.h`）：

```cpp
// 用 lambda 包装旋转谓词
it = std::min_element(first, last,
    [&start_p, &rotation_predicate](const Point& p1, const Point& p2) {
        return rotation_predicate(start_p, p1, p2);
    });

while (!equal_points(*it, stop_p)) {
    *res = *it; ++res;
    it = std::min_element(first, last,
        [it, &rotation_predicate](const Point& p1, const Point& p2) {
            return rotation_predicate(*it, p1, p2);
        });
}
```

**独特接口：`ch_jarvis_march`**

```cpp
// 可以只计算从 start_p 到 stop_p 的部分凸包弧段
ch_jarvis_march(first, last, start_p, stop_p, result, traits);
```

这个接口在其他算法中没有对应物。在某些几何算法中（如计算两个凸包的公切线），只需要凸包的一段弧，Jarvis 是最自然的选择。

**复杂度**：O(n·h)。当 h 很小时优于 O(n log n)；当 h = O(n) 时退化为 O(n²)。

---

### 4.6 Melkman

**角色**：处理简单多边形顶点序列的专用算法，O(n) 时间。

**前提条件**：输入必须是**简单多边形**的顶点序列（按顺序给出，不能是任意点集）。

**核心思路**：用双端队列（`std::deque`）维护当前凸包，在线处理每个新点：

```cpp
std::deque<Point> Q;
// 对每个新点 r：
// 如果 r 在当前凸包内，跳过
// 否则，从队首/队尾弹出不再是凸包顶点的点，插入 r
while (!Q.empty() && !left_turn(r, s, Q.front()))
{ s = Q.front(); Q.pop_front(); }
Q.push_front(s);
```

**复杂度**：O(n)，每个点最多入队出队各一次。

**为什么不能用于任意点集**：Melkman 算法依赖输入点的顺序性——相邻点在多边形边界上相邻。对于任意点集，无法保证这个性质，算法会给出错误结果。

---

## 5. 算法选取的工程逻辑

### 5.1 默认算法的选择标准

`convex_hull_2()` 的选择逻辑：

```
能多次遍历？
  ├─ 是 → ch_akl_toussaint
  │        理由：过滤器效果好，实践中接近线性；底层 Graham-Andrew 可靠
  └─ 否 → ch_bykat
           理由：唯一能处理 input_iterator 的"合理"算法
                 （Eddy 也能，但 Bykat 是改进版）
```

### 5.2 被"放弃"的算法为何仍然实现

| 算法 | 未作为默认的原因 | 保留的原因 |
|------|----------------|-----------|
| Eddy | Bykat 是其改进版，性能更好 | 历史遗留；参考实现；对比测试 |
| Jarvis | O(n·h) 在 h 大时退化 | `ch_jarvis_march` 提供部分凸包弧段接口；h 极小时最优 |
| Melkman | 有严格前提条件，不通用 | O(n) 线性时间；简单多边形场景下无可替代 |
| Graham-Andrew | 被 akl_toussaint 包含 | 作为独立接口暴露；akl_toussaint 的内部构件 |

### 5.3 算法复杂度对比

| 算法 | 时间（最坏） | 时间（平均/实践） | 空间 | 迭代器要求 |
|------|------------|----------------|------|-----------|
| Graham-Andrew | O(n log n) | O(n log n) | O(n) | BidirectionalIterator |
| Akl-Toussaint | O(n log n) | O(n log h) | O(n) | ForwardIterator |
| Bykat | O(n²) | O(n log h) | O(n) | InputIterator |
| Eddy | O(n²) | O(n log h) | O(n) | InputIterator |
| Jarvis | O(n·h) | O(n·h) | O(1) | ForwardIterator |
| Melkman | O(n) | O(n) | O(n) | InputIterator（需简单多边形） |

---

## 6. 后置条件与正确性验证

CGAL 有一套精细的断言宏体系，在 Convex_hull_2 中体现得很典型：

```cpp
// 来自各 _impl.h 文件

// 普通后置条件：验证输出是强凸的（逆时针，无三点共线）
CGAL_postcondition(
    is_ccw_strongly_convex_2(res.output_so_far_begin(),
                             res.output_so_far_end(), ch_traits));

// 昂贵后置条件：暴力验证每个输入点都在凸包内或凸包上
CGAL_expensive_postcondition(
    ch_brute_force_check_2(first, last,
                           res.output_so_far_begin(),
                           res.output_so_far_end(), ch_traits));
```

**`Tee_for_output_iterator`** 是实现后置条件的关键工具：它包装输出迭代器，同时把写出的数据缓存一份，供后置条件检查使用，而不需要修改算法逻辑本身。

```cpp
// 编译期开关：
// CGAL_CH_NO_POSTCONDITIONS → 完全跳过，零开销
// CGAL_NO_POSTCONDITIONS    → 同上
// 默认                      → 启用普通后置条件
// CGAL_EXPENSIVE_POSTCONDITIONS → 启用暴力验证
#if defined(CGAL_CH_NO_POSTCONDITIONS) || defined(CGAL_NO_POSTCONDITIONS)
OutputIterator res(result);
#else
Tee_for_output_iterator<OutputIterator, Point_2> res(result);
#endif
```

这是**零开销抽象**的典型实践：调试时开启验证，发布时编译期消除。

---

## 7. Benchmark 与性能测试

### 7.1 内置 timing 工具

`ch_timing.cpp` 展示了标准的 benchmark 用法：

```cpp
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;

std::vector<Point_2> V(in_start, in_end);
CGAL::ch_timing(V.begin(), V.end(), VE.begin(), iterations, K());
```

`ch_timing` 函数内部对所有算法跑相同的输入，输出各算法的耗时对比。

### 7.2 运行 benchmark

```bash
cd Convex_hull_2/examples/Convex_hull_2
cmake -DCMAKE_BUILD_TYPE=Release .
make ch_timing
./ch_timing files/CD500 100
# 参数：数据文件，迭代次数
```

### 7.3 性能特征与输入分布的关系

| 输入分布 | 最优算法 | 原因 |
|---------|---------|------|
| 均匀随机（圆盘内） | Akl-Toussaint | 过滤器效果极好，h = O(log n) |
| 均匀随机（正方形内） | Akl-Toussaint | 同上 |
| 圆上的点 | Jarvis | h = n，但 Jarvis 的常数小；其他算法 O(n log n) |
| 已排序 | Graham-Andrew | 排序已完成，直接扫描 |
| 流式输入（input_iterator） | Bykat | 唯一选择 |
| 简单多边形顶点 | Melkman | O(n) 无可替代 |

---

## 8. 现代 C++ 模板实践总结

Convex_hull_2 包跨越了 C++98 到 C++14 的演进，代码中可以看到多种模板技术的混用：

### 8.1 Tag Dispatch（C++98 风格，仍然有效）

```cpp
// 通过重载 + tag 对象实现编译期分支
CGAL_convex_hull_points_2(..., std::input_iterator_tag{});
CGAL_convex_hull_points_2(..., std::forward_iterator_tag{});
```

现代替代方案：`if constexpr` + `std::is_base_of_v<std::forward_iterator_tag, Category>`，但 tag dispatch 的可扩展性更好。

### 8.2 Policy-Based Design via Traits

Traits 类是 C++ 泛型编程中的 Policy 模式：算法不依赖具体类型，只依赖 Traits 提供的接口（谓词工厂方法）。这使得同一算法可以在不同精度的 Kernel 下工作，无需修改算法代码。

### 8.3 Lambda 替代函数对象（C++11）

新代码中大量使用 lambda 替代手写 functor：

```cpp
// 旧风格（C++98）：需要单独定义 struct
struct LeftTurnPredicate { ... };

// 新风格（C++11）：就地定义，捕获上下文
auto pred = [&left_turn, &a, &b](const Point_2& p) {
    return left_turn(a, b, p);
};
std::partition(begin, end, pred);
```

### 8.4 std::tuple 返回多值（C++11）

```cpp
// ch_akl_toussaint_impl.h 中返回四个极值点迭代器
std::tuple<ForwardIterator, ForwardIterator, ForwardIterator, ForwardIterator>
ch_nswe_point_with_order(...);

auto [i1, i2, i3, i4] = ch_nswe_point_with_order(...);  // C++17 结构化绑定
```

### 8.5 迭代器特化的层次结构

```
input_iterator_tag
    └── forward_iterator_tag
            └── bidirectional_iterator_tag
                    └── random_access_iterator_tag
```

CGAL 利用这个层次：`forward_iterator_tag` 的重载会匹配所有更强的迭代器类别（因为 C++ 重载解析会选择最匹配的），但 CGAL 选择了显式列出所有类别，避免歧义：

```cpp
// 显式处理每个类别，而不是依赖隐式转换
// 这样代码意图更清晰，也避免了未来新增迭代器类别时的意外行为
```

### 8.6 零开销抽象：编译期可选的后置条件

通过预处理宏在编译期选择是否包含验证代码，发布版本完全没有运行时开销，这是 C++ 相比其他语言的独特优势。

---

## 参考

- CGAL 官方文档：<https://doc.cgal.org/latest/Convex_hull_2/>
- Akl & Toussaint (1978): *A fast convex hull algorithm*
- Andrew (1979): *Another efficient algorithm for convex hulls in two dimensions*
- Bykat (1978): *Convex hull of a finite set of points in two dimensions*
- Jarvis (1973): *On the identification of the convex hull of a finite set of points in the plane*
- Melkman (1987): *On-line construction of the convex hull of a simple polyline*
