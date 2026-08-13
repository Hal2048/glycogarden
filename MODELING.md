# GlycoGarden 糖基化模型建模思路

> 本文档基于 `model_core/` 中的代码，说明 N-糖基化代谢稳态模型的建模思想、数据流、BFS 网络生成机制以及当前状态。完整研究记录位于 Fudan 2026 的独立 Model 仓库。

---

## 1. 问题背景

### 1.1 生物学过程

蛋白质的 **N-连接糖基化（N-glycosylation）** 发生在高尔基体中。初始糖链（通常是 `Man9` 或 `Man8`，即 9 个或 8 个甘露糖）经过高尔基体多个区室时，会被一系列酶逐步修饰：

- **修剪**：ManI、ManII 切除多余甘露糖
- **添加 GlcNAc**：GnTI、GnTII、GnTIII、GnTIV、GnTV、GnTE 等
- **添加半乳糖**：GalT
- **添加岩藻糖**：FucT
- **添加唾液酸**：SiaT

最终形成多种多样的复杂型 N-糖链（glycoforms）。

### 1.2 建模目标

给定：
- 高尔基体各区的停留时间
- 各种酶的总浓度和区室分布
- 糖供体（UDP-GlcNAc、UDP-Gal 等）浓度
- 初始进入的糖链（M9/M8）

预测：
- 稳态时，每种糖链在每个区室的浓度
- 特别是最终区室 TGN 中各种糖链的相对丰度

---

## 2. 模型总体架构

```text
main.py
    │
    ├── config.py          ← 参数：停留时间、酶浓度、供体浓度、酶分布矩阵
    ├── enzyme.py          ← 酶规则：20 条酶反应规则
    ├── reaction.py        ← BFS 生成反应网络
    ├── kinetics.py        ← 竞争抑制计算
    └── golgi_model.py     ← 稳态求解器
                              ↓
                         final_concs[n_structures, 4]
                              ↓
                         取 TGN 列 → 归一化 → top N
```

---

## 3. 各文件功能与联系

### 3.1 `glycoform.py` —— 糖链编码

每种糖链用一个 9 元组表示：

```python
(man, fuc, gnb, br1, br2, br3, br4, gal, sia)
```

| 字段 | 含义 |
|------|------|
| `man` | 甘露糖数量 |
| `fuc` | 岩藻糖数量 |
| `gnb` | 平分型 GlcNAc |
| `br1` ~ `br4` | 4 个分支上的 GlcNAc 状态 |
| `gal` | 半乳糖总数 |
| `sia` | 唾液酸总数 |

每个 `Glycoform` 还维护两个集合：
- `_consuming_reactions`：消耗该糖链的所有反应 ID
- `_producing_reactions`：生成该糖链的所有反应 ID

这两个集合在后续计算净反应速率时非常关键。

### 3.2 `enzyme.py` —— 酶反应规则

定义了 20 条酶规则，每条规则包含：

```python
{
    "name": "酶名",
    "condition": "反应条件表达式",
    "product_struc": "产物结构表达式",
    "cosubstrate": "糖供体",
    "Km": 米氏常数,
    "Kmd": 供体米氏常数,
    "kf": 催化速率常数,
    "adjustments": [("条件", 修正因子), ...]
}
```

`build_enzymes()` 把这些字符串表达式转换成 Python 函数对象。

### 3.3 `reaction.py` —— 反应网络生成（核心：BFS）

这是整个模型最关键的一步。`ReactionNetwork.generate_network()` 使用 **BFS（广度优先搜索）** 从初始糖链生成完整的反应网络。

### 3.4 `kinetics.py` —— 竞争抑制

同一种酶可能作用于多种底物。`KineticsCalculator` 计算每个酶的 competition term：

```
1 + Σ(c_substrate / Km)
```

### 3.5 `golgi_model.py` —— 稳态求解

把高尔基体抽象为 4 个串联 CSTR，建立稳态质量平衡方程，用 `scipy.optimize.root` 求解。

### 3.6 `config.py` —— 参数配置

集中管理所有可调参数，并提供 `build_enzyme_dist(total_enzyme_conc=None, dist_matrix=None)` 把总酶浓度和分布矩阵转换为每个区室的实际酶浓度。

Dashboard 按 10 个基础酶组接收 CGC、MGC、TGC、TGN 四区室的完整分布。`normalize_dist_matrix()` 会验证矩阵形状、有限数值、非负值和正的列总和，并按列归一化，因此每种酶的四区室比例严格合计为 1。`build_enzyme_dist()` 使用请求局部参数，不会为了单次预测修改 `TOTAL_ENZYME_CONC` 或 `DIST_MATRIX` 全局配置。

20 个具体反应酶名称由 10 个基础组展开；例如 GalT 和 SiaT 的各分支反应共享对应基础组的区室分布。这种展开不会改变“每个基础酶组四区室总量等于其总浓度”的守恒约束。

### 3.7 `main.py` —— 执行入口

按顺序调用以上模块，输出 top 糖链。

---

## 4. BFS 反应网络生成详解

### 4.1 为什么用 BFS？

糖链修饰过程天然是一个**状态转移图**：
- 每个节点 = 一种糖链结构
- 每条边 = 一个酶催化反应
- 从初始结构出发，沿着边可以到达很多下游结构

BFS 适合这种**从起点探索所有可达状态**的问题：
- 它能保证不遗漏任何 reachable 的糖链
- 通过队列管理待探索节点
- 配合 `add_structure()` 去重，避免重复处理

### 4.2 代码逐行解析

```python
def generate_network(self, initial_structures: List[Glycoform], enzymes: List[Enzyme]):
    # 1. 把初始糖链加入网络，获取 ID，放入队列
    queue = [self.add_structure(s) for s in initial_structures]

    while queue:
        # 2. 取出队列头部的糖链（FIFO）
        current_id = queue.pop(0)
        current_glycan = self.structures[current_id]

        # 3. 如果这种糖链被标记为 unavailable，跳过
        if current_glycan.available == 0:
            continue

        # 4. 尝试用每一种酶作用于当前糖链
        for enzyme in enzymes:
            if enzyme.condition(current_glycan):
                # 5. 生成产物糖链
                product_glycan = enzyme.product_struc(current_glycan)

                # 6. 把产物加入网络（已存在则返回旧 ID）
                product_id = self.add_structure(product_glycan)

                # 7. 计算这个具体反应的 Km（已调整）
                km_val = enzyme.Km * enzyme.Km_adjust(current_glycan)

                # 8. 创建 Reaction 对象并记录
                rxn_id = self._next_reaction_id
                self._next_reaction_id += 1
                reaction = Reaction(
                    rxn_id=rxn_id,
                    enzyme=enzyme,
                    sub_id=current_id,
                    prod_id=product_id,
                    km=km_val
                )
                self.reactions[rxn_id] = reaction

                # 9. 更新糖链的消耗/生成关系
                current_glycan._consuming_reactions.add(rxn_id)
                product_glycan._producing_reactions.add(rxn_id)
                enzyme.substrate_km_map[current_id] = km_val

                # 10. 如果产物还没被探索过，加入队列
                if product_id not in queue:
                    queue.append(product_id)
```

### 4.3 BFS 执行示例

假设初始糖链是 `Man9`，酶规则简化为：

- ManI：如果 `man > 5`，产物 `man-1`
- ManII：如果 `man > 3 and br4 == 1`，产物 `man-1`

执行过程：

```text
初始队列: [Man9]

第 1 轮:
  取出 Man9
  ManI 作用 → Man8, 队列: [Man8]

第 2 轮:
  取出 Man8
  ManI 作用 → Man7, 队列: [Man7]

第 3 轮:
  取出 Man7
  ManI 作用 → Man6, 队列: [Man6]

第 4 轮:
  取出 Man6
  ManI 作用 → Man5, 队列: [Man5]

第 5 轮:
  取出 Man5
  ManI 不作用（man 不大于 5）
  假设 GnTI 作用，生成 br4=1 的结构，队列: [新结构]

...继续直到队列为空
```

### 4.4 关键设计点

#### 去重机制

```python
def add_structure(self, glycan: Glycoform) -> int:
    if glycan in self._structure_to_id:
        return self._structure_to_id[glycan]
    struct_id = self._next_id
    self._next_id += 1
    self.structures[struct_id] = glycan
    self._structure_to_id[glycan] = struct_id
    return struct_id
```

同一种产物可能由不同底物生成，但只存储一次。返回的是唯一 ID。

#### 为什么用集合存消耗/生成反应？

```python
current_glycan._consuming_reactions.add(rxn_id)
product_glycan._producing_reactions.add(rxn_id)
```

在 `golgi_model.py` 中计算净速率时：

```python
for i, glycan in self.network.structures.items():
    for rxn_id in glycan._consuming_reactions:
        r[i] -= rxn_rates.get(rxn_id, 0.0)
    for rxn_id in glycan._producing_reactions:
        r[i] += rxn_rates.get(rxn_id, 0.0)
```

这样可以直接根据糖链找到所有相关反应，O(1) 访问。

#### 终止条件

队列为空时停止。什么情况下会空？
- 所有 reachable 糖链都被处理过
- 新产物都是已存在的结构，不再加入队列
- 或者产物 `available == 0` 被跳过

#### 规模

实际运行中：
- 初始：2 种（Man9、Man8）
- 最终：约 **7,500 种结构**
- 反应：约 **22,800 个**

### 4.5 BFS 的局限性

- 规则必须保证网络不会无限增长（通过 `available` 限制分支数）
- 如果酶规则有循环（A→B→A），可能永远跑不完。当前规则设计是单向修饰，不会出现循环
- 内存消耗随网络规模线性增长

---

## 5. 稳态求解详解

### 5.1 模型假设

把高尔基体 4 个区室看成 **串联的 CSTR（连续搅拌反应器）**：

```
Feed → CGC → MGC → TGC → TGN
```

每个区室中，糖链浓度达到稳态：

```
c_{i,j} = c_{i,j-1} + τ_j * r_{i,j}
```

- `c_{i,j}`：糖链 `i` 在区室 `j` 的浓度
- `τ_j`：在区室 `j` 的停留时间
- `r_{i,j}`：糖链 `i` 在区室 `j` 的净生成速率

### 5.2 反应速率

单个反应速率用修正的 Michaelis-Menten 方程：

```
r = kf * [Et] * [donor] * [Pi] / (Kmi * (Kmd + [donor]) * (1 + Σ(Pj)/Kmj))
```

### 5.3 求解过程

1. 初始化猜测值 `x0`：每个区室都等于初始 feed
2. 计算残差 `F(x) = c - (c_prev + τ * r)`
3. 用 `scipy.optimize.root` 找到 `F(x) = 0`
4. 返回 `final_concs`

### 5.4 输出

```python
final_concs.shape == (n_structures, 4)
```

最后一列 `final_concs[:, -1]` 就是 TGN 区室的最终产物分布。

---

## 6. 输入与输出

### 6.1 输入

| 输入 | 文件 | 含义 |
|------|------|------|
| `TAU` | config.py | 4 个区室停留时间 |
| `DONOR_CONC` | config.py | 糖供体浓度 |
| `TOTAL_ENZYME_CONC` | config.py | 10 类酶总浓度 |
| `DIST_MATRIX` | config.py | 酶在 4 个区室的分布比例 |
| 初始 feed | main.py | M9 和 M8 各 50% |

### 6.2 输出

- `final_concs`：每种糖链在 4 个区室的稳态浓度
- TGN 列：最终产物分布
- Top N 糖链：按浓度排序的前 N 种结构

---

## 7. Software 中的应用

`backend/api.py` 将前端的具名分布转换为模型矩阵：

```json
{
  "ManI": {"CGC": 0.05, "MGC": 0.15, "TGC": 0.40, "TGN": 0.40}
}
```

实际请求必须包含配置接口返回的全部 10 个基础酶组。后端拒绝缺失或未知的酶/区室、布尔值、非数字、非有限数、负数和四项全零的分布，然后再次归一化并将实际采用的 `enzymeDistribution` 放入预测响应。

用户界面（`frontend/`）提供：

- 启动子强度选择（4 级）；
- CGC、MGC、TGC、TGN 四区室直接输入；
- 按比例联动、区室锁定、单酶重置和全部重置；
- Baseline、CGC-biased、TGN-biased、Uniform 完整矩阵预设；
- 对应条件下的 top 糖链柱状图和表格。

任意手动分布必须调用实时模型，不能由有限的预计算组合覆盖。`generate_matrix()` 仅生成启动子 × 预设的离线结果，当前格式版本为 `2.0.0`。

---

## 8. 当前状态与计算限制

### 8.1 已完成的

- 糖链编码方案确定
- 20 条酶规则定义
- BFS 反应网络生成
- 竞争抑制计算
- 稳态求解器框架
- Software 四区室交互与严格分布校验
- 请求局部的酶浓度/分布参数，不修改模型全局配置

### 8.2 计算限制

完整网络约有 7,500 种糖链结构，单次稳态求解可能需要约 1–2 分钟。因此页面在计算期间禁止重复提交，并在参数改变后将旧结果标记为过期。单个后端进程会串行执行求解，并缓存最近 32 个参数完全相同的已完成结果。页面的“停止等待”只中止浏览器请求，不会伪装成服务器端求解取消；服务器仍可能完成计算并写入缓存。分布式任务队列和跨进程缓存属于后续运行架构优化。

旧版 `matrix.json` 使用 `1.0.0` schema 和 Early/Late bias，不能表达新的完整四区室预设。需要显式运行 `python backend/generate_matrix.py` 生成 `2.0.0` 数据；服务不会静默返回旧格式结果。

---

## 9. 总结

这个模型的核心思路是：

1. **编码**：用 9 元组数字化每种糖链
2. **规则**：用 20 条酶规则定义糖链如何转化
3. **BFS 建网**：从 M9/M8 出发，遍历所有 reachable 糖链和反应
4. **动力学**：用修正 Michaelis-Menten 计算每个反应速率
5. **稳态求解**：解 4-CSTR 质量平衡方程
6. **输出**：得到最终糖链分布

BFS 是整个流程的枢纽：它把离散的酶规则和初始结构，扩展成完整的反应网络，为后续的动力学和稳态计算提供基础。
