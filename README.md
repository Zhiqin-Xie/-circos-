# 基因组相似性circos弦图绘制工具
物种基因组丰度相似性对比工具，可以将需要对比的基因组文件（.fna格式）和对比目标基因组文件分别放置在不同文件夹，通过运行脚本，生成对比图（png，svg，PDF格式）以及ANI对比表格（csv格式）。另有两个json格式文件，记录我们的基因组之间相似度大小以及相似位置长短以及偏移位置；生成json的另一个目的是可以在脚本中修改简单的绘图参数后，通过json文件快速重绘图片，避免微小修改耗费大量时间重绘。

# 基因组Circos可视化脚本 - 技术原理与工作流程

> 基因组Circos可视化-从环境配置到图形生成的完整技术解析

---

## 📋 目录

1. [环境配置](#环境配置)
2. [文件准备](#文件准备)
3. [脚本运行流程](#脚本运行流程)
4. [计算原理](#计算原理)
5. [可视化原理](#可视化原理)
6. [缓存机制](#缓存机制)

---

## 1. 环境配置

### 1.1 必需工具

#### Python环境（conda推荐）

```bash
# 创建conda环境
conda create -n circos python=3.8

# 激活环境
conda activate circos
```

#### Python依赖包

```bash
# 核心可视化包
pip install pycirclize

# 数据处理包
pip install pandas numpy matplotlib

# 或使用conda安装
conda install pandas numpy matplotlib
conda install -c conda-forge pycirclize
```

#### 生物信息学工具

```bash
# 必需：MUMmer（包含nucmer）
conda install -c bioconda mummer

# 可选但推荐：fastANI
conda install -c bioconda fastani
```

### 1.2 环境验证

```bash
# 验证工具是否安装成功
nucmer --version     # 应显示版本号
fastANI -h          # 应显示帮助信息
python -c "import pycirclize; print('OK')"
```

---

## 2. 文件准备

### 2.1 文件夹结构

```
工作目录/
├── genome_circos_qva_fast.py    # 脚本文件
├── reference/                    # 参考基因组文件夹
│   ├── genome1.fna
│   ├── genome2.fna
│   └── ...
├── query/                        # 查询基因组文件夹
│   ├── target1.fna
│   ├── target2.fna
│   └── ...
└── run05/                        # 输出目录（自动创建）
    ├── circos_comparison.png
    ├── ani_matrix.csv
    └── *_cache.json
```

### 2.2 创建文件夹

```bash
# 进入工作目录
cd /path/to/your/work/directory

# 创建必需文件夹
mkdir -p reference query

# 确认结构
tree -L 1
# 应显示：
# .
# ├── genome_circos_qva_fast.py
# ├── reference/
# └── query/
```

### 2.3 放置基因组文件

#### 文件格式要求

- **格式**：FASTA格式（.fna, .fasta, .fa）
- **内容**：组装完成的基因组序列
- **完整性**：可以是完整基因组或草图

#### 放置规则

```bash
# Reference文件夹：普通/背景基因组
cp /path/to/normal/genomes/*.fna reference/

# Query文件夹：重点关注的基因组
cp /path/to/target/genomes/*.fna query/

# 示例：
# reference/ - 45个常见菌株
# query/     - 5个新分离的菌株
```

#### 命名建议

```
推荐格式：GCA_000123456.1_organism_name.fna
        └─ 登录号   └─ 描述性名称

或简化格式：strain01.fna, strain02.fna

脚本会自动简化名称：
  GCA_000508745.1_Stp_DORA_6_22_velvet.fna
  → 显示为：GCA_000508745.1
```

---

## 3. 脚本运行流程

### 3.1 完整运行流程（mode=full）

```
启动脚本
   ↓
1. 检查工具可用性
   ├─ nucmer（必需）
   └─ fastANI（可选）
   ↓
2. 读取基因组文件
   ├─ reference/*.fna
   └─ query/*.fna
   ↓
3. 计算ANI矩阵（如果有fastANI）
   ├─ 每对基因组计算一次
   ├─ 结果：相似度百分比
   └─ 保存：ani_matrix.csv
   ↓
4. 执行全基因组比对（nucmer）
   ├─ Query vs All策略
   ├─ 每对基因组执行一次nucmer
   ├─ 解析相似区域
   └─ 过滤：按相似度和长度
   ↓
5. 保存缓存文件
   ├─ links_cache.json（连线数据）
   └─ genomes_info.json（基因组信息）
   ↓
6. 绘制Circos图
   ├─ 构建圆形布局
   ├─ 绘制基因组块
   ├─ 绘制标签
   ├─ 绘制连线
   └─ 保存图片（PNG/PDF/SVG）
   ↓
完成
```

### 3.2 快速重绘流程（mode=plot）

```
启动脚本
   ↓
1. 检查缓存文件
   ├─ links_cache.json ✓
   └─ genomes_info.json ✓
   ↓
2. 加载缓存数据
   ├─ 跳过ANI计算
   ├─ 跳过nucmer比对
   └─ 直接读取连线数据
   ↓
3. 绘制Circos图（几十秒）
   ├─ 应用新的视觉参数
   ├─ 重新渲染
   └─ 保存图片
   ↓
完成
```

---

## 4. 计算原理

### 4.1 fastANI：平均核苷酸一致性

#### 原理

fastANI使用**基于k-mer的快速比对算法**计算全基因组平均核苷酸一致性（Average Nucleotide Identity, ANI）。

```
基因组A ──┐
          ├─→ [k-mer匹配] ─→ [加权平均] ─→ ANI值（%）
基因组B ──┘
```

#### 计算过程

```
步骤1: 基因组分片
  基因组 → 1kb片段集合

步骤2: 快速映射
  片段A → 寻找基因组B中的最佳匹配位置

步骤3: 计算相似度
  每个片段 → 计算与匹配区域的序列一致性

步骤4: 加权平均
  所有片段的一致性 → 加权平均 → ANI值
```

#### 输出结果

```
基因组A vs 基因组B
  ANI = 98.5%    # 平均相似度

解读：
  > 95%    同一种（species）
  90-95%   可能同一种
  75-90%   同一属（genus）
  < 75%    计算失败（差异太大）
```

#### 在脚本中的作用

```python
# 生成相似度矩阵（N×N）
ani_matrix = calculate_all_ani(files)

# 结果示例：
#              GCA_001  GCA_002  GCA_003
# GCA_001       100.0     98.5     76.3
# GCA_002        98.5    100.0     75.8
# GCA_003        76.3     75.8    100.0
```

**注意**：

- ANI计算是**全局相似度**
- 结果用于参考，不直接影响Circos图的连线
- 连线是由nucmer的**局部比对**决定的

---

### 4.2 nucmer：核苷酸精确比对

#### 原理

nucmer（NUCleotide MUMmer）使用**后缀树算法**进行全基因组精确比对，找出所有相似区域。

```
基因组A: ════════════════════════════
基因组B: ════════════════════════════
         ↓ nucmer比对
找到相似区域：
  区域1: A[100k-150k] ←→ B[200k-250k], 98%相似
  区域2: A[300k-350k] ←→ B[450k-500k], 92%相似
  区域3: A[800k-850k] ←→ B[100k-150k], 88%相似
```

#### 比对步骤

```
步骤1: 构建后缀树
  基因组A → 后缀树索引

步骤2: 查找最大唯一匹配（MUM）
  基因组B → 在A的后缀树中查找匹配片段

步骤3: 扩展匹配区域
  MUM → 向两端扩展 → 完整比对区域

步骤4: 计算比对参数
  每个区域：
    - 起始位置（ref_start, query_start）
    - 结束位置（ref_end, query_end）
    - 序列一致性（identity %）
    - 区域长度（length bp）
    - 错配数（errors）
```

#### 输出：Delta文件

```
>genome_A_scaffold1 >genome_B_scaffold1 2814816 2750948
100000 150000 200000 250000 500 0 0
300000 350000 450000 500000 800 0 0
...

解读：
  第一行：scaffold信息和长度
  后续行：每个比对区域
    - 前4个数字：起止坐标
    - 第5个数字：错配数
```

#### 过滤参数

脚本中的两个关键参数：

```python
MIN_IDENTITY = 85.0   # 相似度阈值
MIN_LENGTH = 500      # 长度阈值

# 过滤逻辑
for region in all_regions:
    identity = (1 - errors/length) * 100
    if identity >= MIN_IDENTITY and length >= MIN_LENGTH:
        keep_region()  # 保留用于绘制连线
    else:
        discard_region()  # 丢弃
```

#### 结果数据结构

```python
links = [
    {
        'ref_start': 100000,      # 基因组A起始位置
        'ref_end': 150000,        # 基因组A结束位置
        'query_start': 200000,    # 基因组B起始位置
        'query_end': 250000,      # 基因组B结束位置
        'identity': 98.5,         # 相似度百分比
        'length': 50000           # 区域长度（bp）
    },
    ...
]
```

---

### 4.3 两种计算的关系

```
fastANI（全局）:
  计算整个基因组的平均相似度
  用途：快速评估总体相似性
  输出：单个百分比值
  速度：快（10-30秒/对）

nucmer（局部）:
  找出所有相似区域的具体位置
  用途：精确定位相似片段
  输出：多个区域的详细信息
  速度：较慢（5-15秒/对）

关系：
  fastANI告诉你"总体多相似"
  nucmer告诉你"哪里相似"

  ANI = 98% ↗
              → 基因组很相似
  nucmer找到100个区域 ↗

  ANI计算失败（<75%）↗
              → 基因组差异大
  nucmer可能只找到几个区域 ↗
```

---

## 5. 可视化原理

### 5.1 Circos图基本原理

#### 坐标系统

```
径向坐标（r）：
  0（圆心）→ 100（外圈）

  标签：r = 99-130（圆外）
  基因组块：r = 92-98（外圈）
  连线：r = 0-92（内部）

角度坐标（θ）：
  0° - 360°（顺时针）
  每个基因组：起始角 → 结束角
  间隔：SECTOR_SPACE（度）
```

### 5.2 绘制流程

#### 步骤1：构建Sectors（基因组块）

```python
# 计算每个基因组的弧长
sectors = {
    'genome1': length1,  # 2,814,816 bp
    'genome2': length2,  # 2,750,948 bp
    ...
}

# 转换为角度
total_length = sum(all_lengths)
total_angle = 360 - (n_genomes × SECTOR_SPACE)

genome1_angle = (length1 / total_length) × total_angle
genome2_angle = (length2 / total_length) × total_angle
```

#### 步骤2：绘制基因组块

```python
for genome in genomes:
    # 外圈彩色轨道
    track = add_track(r_inner=92, r_outer=98)
    track.axis(color=genome.color)

    # 颜色编码
    if genome.is_query:
        color = warm_colors  # 红橙黄
        thickness = 2.5
    else:
        color = cold_colors  # 蓝绿青
        thickness = 1.5
```

#### 步骤3：绘制标签

```python
# 径向标签（LABEL_ADJUST_ROTATION=True）
label.position = (angle, r=125)
label.rotation = angle - 90°  # 切向排列

# 固定方向标签（LABEL_ADJUST_ROTATION=False）
label.position = (angle, r=110)
label.rotation = 90°  # 全部竖直
```

#### 步骤4：绘制连线（Bezier曲线）

```python
# 每个相似区域 → 一条连线
for link in all_links:
    # 起点：基因组A的区域
    start = (genomeA, start_pos, end_pos)

    # 终点：基因组B的区域
    end = (genomeB, start_pos, end_pos)

    # 绘制贝塞尔曲线
    bezier_curve(start, end, 
                 color=genomeA.color,
                 alpha=f(identity),  # 透明度随相似度变化
                 width=f(length))    # 粗细随长度变化
```

### 5.3 视觉编码

#### 颜色

```
Reference基因组：冷色调
  蓝色 → 青色 → 绿色
  HSV: Hue=180°-270°

Query基因组：暖色调
  红色 → 橙色 → 黄色
  HSV: Hue=0°-60°
```

#### 连线属性

```
透明度（Alpha）：
  基础值：query=0.4, ref=0.25
  相似度加成：(identity - MIN_IDENTITY) / 50 × 0.2
  最大值：0.6

  效果：相似度越高，连线越不透明

粗细（Width）：
  query连线：0.9
  reference连线：0.6

  效果：query的连线更突出
```

#### 基因组块

```
厚度：
  track_outer - track_inner

  query：98 - 92 = 6个单位
  ref：  98 - 95 = 3个单位

  效果：query基因组块更厚（更突出）

边框：
  query：linewidth=2.5（粗边框）
  ref：  linewidth=1.5（细边框）
```

---

## 6. 缓存机制

### 6.1 缓存文件

#### links_cache.json（连线数据）

```json
[
  {
    "g1_name": "GCA_000508745.1",
    "g2_name": "GCA_001914135.1",
    "g1_is_ref": false,
    "links": [
      {
        "ref_start": 100000,
        "ref_end": 150000,
        "query_start": 200000,
        "query_end": 250000,
        "identity": 98.5,
        "length": 50000
      },
      ...
    ]
  },
  ...
]
```

#### genomes_info.json（基因组信息）

```json
[
  {
    "name": "GCA_000508745.1",
    "short_name": "GCA_000508745.1",
    "file": "reference/GCA_000508745.1.fna",
    "length": 2814816,
    "color": "#4a90e2",
    "is_reference": true
  },
  ...
]
```

### 6.2 工作原理

```
首次运行（mode=full）:
  1. 执行所有计算
  2. 生成结果数据
  3. 保存到缓存文件
  4. 绘制图形

后续运行（mode=plot）:
  1. 检测到缓存文件存在
  2. 跳过ANI计算（节省2-4小时）
  3. 跳过nucmer比对（节省1-2小时）
  4. 直接从缓存读取数据
  5. 应用新的视觉参数
  6. 重新绘制图形（30秒）
```

### 6.3 缓存失效条件

```
需要重新运行mode=full的情况：

1. 数据改变
   - 添加/删除基因组文件
   - 更换基因组序列

2. 过滤参数改变
   - MIN_IDENTITY变化
   - MIN_LENGTH变化

3. 缓存文件丢失
   - links_cache.json被删除
   - genomes_info.json被删除
```

---

## 7. 完整工作示例

### 7.1 实际运行示例

#### 环境准备

```bash
# 1. 激活环境
conda activate circos

# 2. 进入工作目录
cd ~/genome_circos_test

# 3. 确认文件结构
tree -L 2
# .
# ├── genome_circos_qva_fast.py
# ├── reference/
# │   ├── genome1.fna
# │   └── ...（45个文件）
# └── query/
#     ├── target1.fna
#     └── ...（5个文件）
```

#### 首次运行

```bash
python genome_circos_qva_fast.py --mode full > run.log 2>&1 &

# 实时监控进度
tail -f run.log

# 输出示例：
# Checking tools...
#   nucmer: ✓
#   fastANI: ✓
# 
# 计算ANI值
# [1/1225] GCA_001 vs GCA_002
#          ANI = 98.5%
# ...
# 
# 执行全基因组比对 (Query vs All)
# [1/245] [Query] target1 <-> [Ref] genome1
#   -> Found 15 similar regions
# ...
# 
# 绘制Circos图
# Total links drawn: 850
# Saving images...
# ✓ All tasks completed!
```

#### 计算过程时间分配

```
50个基因组（5 query + 45 ref）示例：

ANI计算：
  1,225对 × 15秒 = 5.1小时

nucmer比对：
  245对 × 10秒 = 40分钟

绘图：
  数十秒

总计：约6小时
```

#### 快速重绘

```bash
# 修改参数
nano genome_circos_qva_fast.py
# 改：LABEL_RADIUS = 125

# 快速重绘
python genome_circos_qva_fast.py --mode plot

# 输出：
# Loading cache files...
#   ✓ Loaded 245 link groups
#   ✓ Loaded 50 genomes
# 
# 绘制Circos图
# Total links drawn: 850
# Saving images...
# ✓ Circos plot completed!
```

### 7.2 数据流转图

```
输入文件
  ├─ reference/*.fna (45个)
  └─ query/*.fna (5个)
        ↓
   读取基因组
        ↓
     计算ANI ──────────────→ ani_matrix.csv
        ↓
    nucmer比对
        ↓
   解析delta文件
        ↓
    过滤区域
   (identity≥85%, length≥500)
        ↓
    保存缓存 ──────────────→ links_cache.json
        ↓                    genomes_info.json
  构建Circos布局
        ↓
   绘制基因组块
   （颜色、厚度、边框）
        ↓
    绘制标签
   （位置、角度、字体）
        ↓
    绘制连线
  （贝塞尔曲线、透明度）
        ↓
    渲染图形 ──────────────→ circos_comparison.png
        ↓                    circos_comparison.pdf
    保存文件                 circos_comparison.svg
```

---

## 8. 技术要点总结

### 8.1 核心算法

| 工具  | 算法  | 输出  | 用途  |
| --- | --- | --- | --- |
| fastANI | k-mer匹配 | 全局相似度(%) | 快速评估 |
| nucmer | 后缀树+MUM | 局部相似区域 | 精确定位 |
| pycirclize | 极坐标绘图 | Circos图 | 可视化 |

### 8.2 数据规模

```
50个基因组（每个~3Mb）:

ANI计算：
  - 配对数：50×49/2 = 1,225
  - 内存：~100MB
  - 时间：~5小时

nucmer比对：
  - 配对数：5×49 = 245（query vs all）
  - 内存：~500MB/对
  - 时间：~40分钟
  - 临时文件：~50GB（会自动清理）

最终图形：
  - 内存：~2GB
  - 输出：~10MB（PNG）
```

### 8.3 性能优化

```
1. 缓存机制
   - 避免重复计算
   - 几十秒快速重绘

2. 并行计算
   - nucmer自动多线程
   - ANI可串行（已够快）

3. 内存管理
   - 逐对计算不积累
   - 及时清理临时文件
```

### 8.4 关键参数影响

| 参数  | 增大效果 | 减小效果 |
| --- | --- | --- |
| MIN_IDENTITY | 连线更少（更严格） | 连线更多（更宽松） |
| MIN_LENGTH | 只保留长片段 | 包含短片段 |
| SECTOR_SPACE | 基因组块更分散 | 基因组块更紧凑 |
| LABEL_RADIUS | 标签更远离圆 | 标签更靠近圆 |
| LINK_ALPHA | 连线更不透明 | 连线更透明 |

---

## 9. 常见错误故障排查

#### 错误1：nucmer not found

```bash
# 原因：未安装MUMmer
# 解决：
conda install -c bioconda mummer
```

#### 错误2：plot mode requires cache

```bash
# 原因：没有缓存文件
# 解决：
python genome_circos_qva_fast.py --mode full
```

---

## 10. 总结

### 工作流程精简版

```
准备 → 计算 → 缓存 → 绘图

详细：
1. 配置环境（conda + 工具）
2. 准备基因组文件（reference + query）
3. 运行mode=full（3-6小时）
   - ANI计算（全局相似度）
   - nucmer比对（局部区域）
4. 保存缓存（JSON文件）
5. 绘制Circos图（数十秒）
6. 后续执行脚本时用命令调整参数用mode=plot重绘（数十秒）
```
