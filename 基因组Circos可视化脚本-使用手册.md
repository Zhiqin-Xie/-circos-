# 基因组Circos可视化工具 - 使用手册

> Query vs All 比对策略 + 快速重绘版本

---

## 一、核心功能

- **Query vs All 比对**：将数据分别放置于两个文件夹（query和reference），只有query文件夹中的基因组数据会与所有基因组（query和reference）比对
- **快速重绘**：调整视觉参数后数十秒出图（无需重新计算）
- **缓存机制**：自动保存比对结果，支持多次调整
- **高质量输出**：PNG/PDF/SVG三种格式

---

## 二、快速开始

### 1. 准备数据

```bash
# 创建文件夹
mkdir -p reference query

# 放入基因组文件
cp normal_genomes/*.fna reference/    # 普通基因组、参考基因组
cp target_genomes/*.fna query/        # 重点基因组（会突出显示）
```

### 2. 第一次运行

```bash
python genome_circos_qva_fast.py --mode full > run.log 2>&1 &

# 等待3-6小时...
# 生成图片 + 缓存文件（用于部分绘图相关参数修改后使用plot模式快速重绘，后续有详细说明）
```

### 3. 调整参数重绘

```bash
# 编辑脚本修改参数
nano genome_circos_qva_fast.py

# 快速重绘
python genome_circos_qva_fast.py --mode plot > run.log 2>&1 &
```

---

## 三、参数配置

### 1.需要重新计算的参数（mode=full）

修改这些参数后必须用 `--mode full` 重新运行：

| 参数  | 说明  | 默认值 |
| --- | --- | --- |
| `REFERENCE_DIR` | 参考基因组目录 | `"reference"` |
| `QUERY_DIR` | 查询基因组目录 | `"query"` |
| `NUM_QUERY_GENOMES` | 选择几个query基因组 | `-1`（-1默认将query文件夹中所有基因组用于和reference文件夹中的基因组对比） |
| `MIN_IDENTITY` | 相似度过滤阈值（%） | `85.0` |
| `MIN_LENGTH` | 片段长度过滤阈值（bp） | `500` |

**改这些参数后的操作（需要重新进行ANI计算和）**：

```bash
python genome_circos_qva_fast.py --mode full
```

---

### 2.只需重绘的参数（mode=plot）

修改这些参数后只需 `--mode plot`，快速出图：

#### 布局参数

| 参数  | 说明  | 推荐值 |
| --- | --- | --- |
| `SECTOR_SPACE` | 基因组块间隔（度） | `0.5-2` |
| `LABEL_RADIUS` | 标签半径 | `99-130` |
| `GENOME_TRACK_INNER_QUERY` | query基因组块内径 | `92` |
| `GENOME_TRACK_INNER_REF` | reference基因组块内径 | `93-95` |
| `GENOME_TRACK_OUTER` | 基因组块外径 | `98` |

#### 标签样式

| 参数  | 说明  | 推荐值 |
| --- | --- | --- |
| `LABEL_SIZE_QUERY` | query标签字体大小 | `6-10` |
| `LABEL_SIZE_REF` | reference标签字体大小 | `5-8` |
| `LABEL_ADJUST_ROTATION` | 标签是否径向排列 | `True` |
| `LABEL_ORIENTATION` | 标签方向 | `None`  见下方说明 |

#### 连线样式

| 参数  | 说明  | 推荐值 |
| --- | --- | --- |
| `LINK_ALPHA_BASE_QUERY` | query连线透明度 | `0.3-0.6` |
| `LINK_ALPHA_BASE_REF` | reference连线透明度 | `0.2-0.5` |
| `LINK_WIDTH_QUERY` | query连线粗细 | `0.6-1.5` |
| `LINK_WIDTH_REF` | reference连线粗细 | `0.4-1.0` |

#### 图片输出

| 参数  | 说明  | 推荐值 |
| --- | --- | --- |
| `FIGURE_SIZE` | 图片大小（英寸） | `(18,18)-(30,30)` |
| `FIGURE_DPI` | 分辨率 | `400-1200` |

**改这些参数后的操作，可快速重新绘图，避免了耗时的ANI计算和nucmer比对**：

```bash
python genome_circos_qva_fast.py --mode plot
```

---

## 四、常见调整场景

### 场景1：标签、物种名重叠

**问题**：基因组多时标签挤在一起

**解决方案**：设置径向标签

```python
# 找到第70-71行，修改为：
LABEL_ADJUST_ROTATION = True   # 启用径向排列
LABEL_ORIENTATION = None       # 注意：是 None 不是 "None"
```

**重要**：`LABEL_ORIENTATION` 必须是 `None`（无引号），不是字符串 `"None"`

---

### 场景2：基因组块太小

**解决方案**：减小间隔

```python
SECTOR_SPACE = 0.5    # 从 0.8度 改为 0.5度
```

```bash
python genome_circos_qva_fast.py --mode plot
```

---

### 场景3：连线太淡看不清

**解决方案**：增加透明度

```python
LINK_ALPHA_BASE_QUERY = 0.6   # 从 0.4 改为 0.6
LINK_ALPHA_BASE_REF = 0.4     # 从 0.25 改为 0.4
```

```bash
python genome_circos_qva_fast.py --mode plot
```

---

### 场景4：需要更高分辨率

**解决方案**：增大尺寸和DPI

```python
FIGURE_SIZE = (30, 30)   # 从 (23,23) 改为 (30,30)
FIGURE_DPI = 900        # 从 700 改为 900
```

```bash
python genome_circos_qva_fast.py --mode plot
```

---

### 场景5：想要更多/更少连线

**解决方案**：调整过滤阈值（需要重新计算）

```python
# 更多连线（宽松）
MIN_IDENTITY = 75.0    # 从 85.0 改为 75.0
MIN_LENGTH = 500       # 保持不变

# 更少连线（严格）
MIN_IDENTITY = 90.0    # 从 85.0 改为 90.0
MIN_LENGTH = 10000     # 从 500 改为 10000
```

```bash
python genome_circos_qva_fast.py --mode full  # 需要 full 模式
```

---

## 五、运行模式

### 模式1：full（完整运行）

```bash
python genome_circos_qva_fast.py --mode full

# 用途：第一次运行，或改变数据/过滤参数时
# 耗时：3-6小时
# 生成：图片 + ANI矩阵 + 缓存文件
```

---

### 模式2：plot（快速重绘）

```bash
python genome_circos_qva_fast.py --mode plot

# 用途：调整视觉参数后快速重绘
# 耗时：30秒！
# 前提：必须先运行过 full 模式
```

---

### 模式3：auto（自动检测，推荐）

```bash
python genome_circos_qva_fast.py

# 自动判断：
#   - 有缓存文件 → 使用 plot 模式（30秒）
#   - 无缓存文件 → 使用 full 模式（3-6小时）
```

---

## 六、输出文件

```
run/
├─ circos_comparison.png    # 高分辨率PNG图
├─ circos_comparison.pdf    # 矢量PDF图
├─ circos_comparison.svg    # 可编辑SVG图
├─ ani_matrix.csv           # ANI相似度矩阵
├─ links_cache.json         # 连线数据缓存
└─ genomes_info.json        # 基因组信息缓存
```

**重要**：不要删除 `*_cache.json` 文件，它们是快速重绘的关键！

---

## 七、常见问题

### Q1: 标签重叠怎么办？

**A**: 设置径向标签

```python
LABEL_ADJUST_ROTATION = True
LABEL_ORIENTATION = None    # 不是 "None"
```

然后运行 `--mode plot`

---

### Q2: 改了参数后用哪个模式？

**A**: 看参数类型

- 改了 **数据源** 或 **过滤参数** → `--mode full`（3-6小时）
- 改了 **布局/样式/输出** 参数 → `--mode plot`（30秒）

---

### Q3: 提示 "plot mode requires cache files"

**A**: 缺少缓存文件，需要先运行一次：

```bash
python genome_circos_qva_fast.py --mode full
```

---

### Q4: 如何强制重新计算？

**A**: 删除缓存或使用 full 模式

```bash
# 方法1：删除缓存
rm run05/*_cache.json

# 方法2：直接用 full
python genome_circos_qva_fast.py --mode full
```

---

## 八、首次运行

```bash
# 1. 准备数据
mkdir -p reference query
cp *.fna reference/ #参考基因组
cp *.fna query/ #需要比对的基因组

# 2. 完整运行
python genome_circos_qva_fast.py --mode full > run.log 2>&1 &

# 3. 监控进度
tail -f run.log
```

---

## 九、参数速查表

### 一张表看懂所有参数

| 参数类别 | 参数名 | 改后使用模式 | 常用值 |
| --- | --- | --- | --- |
| 数据  | `REFERENCE_DIR`, `QUERY_DIR` | **full** | 文件夹路径 |
| 数据  | `NUM_QUERY_GENOMES` | **full** | `-1` 或数字 |
| 过滤  | `MIN_IDENTITY` | **full** | `75-95` |
| 过滤  | `MIN_LENGTH` | **full** | `500-50000` |
| 布局  | `SECTOR_SPACE` | **plot** | `0.5-2` |
| 布局  | `LABEL_RADIUS` | **plot** | `99-130` |
| 布局  | `GENOME_TRACK_*` | **plot** | `60-98` |
| 标签  | `LABEL_SIZE_*` | **plot** | `5-10` |
| 标签  | `LABEL_ADJUST_ROTATION` | **plot** | `True/False` |
| 标签  | `LABEL_ORIENTATION` | **plot** | `None/vertical` |
| 连线  | `LINK_ALPHA_*` | **plot** | `0.2-0.8` |
| 连线  | `LINK_WIDTH_*` | **plot** | `0.4-2.0` |
| 输出  | `FIGURE_SIZE` | **plot** | `(18,18)-(30,30)` |
| 输出  | `FIGURE_DPI` | **plot** | `400-1200` |

**plot = 快速重绘**

---

## 十、核心要点

1. **两类参数**：计算参数（需full）vs 视觉参数（只需plot）
2. **快速迭代**：视觉参数随便改，可快速重绘见效果
3. **保护缓存**：不要删除 `*_cache.json`
4. **auto模式**：最省心，自动判断
