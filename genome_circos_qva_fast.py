#!/usr/bin/env python3
"""
基因组相似性Circos图可视化工具 - Query vs All

核心特性：
===============================================
1. Query vs All 比对策略
2. 快速重绘机制 (支持 --mode plot 快速出图)
3. 中间结果缓存 (记录包含全局Offset的精准坐标)
4. 多Contig基因组支持 (完美解决连线堆叠问题)
5. 自动清理 nucmer 残留的临时文件

使用方法：
===============================================
# 第一次运行（完整计算）
python genome_circos_qva_fast.py --mode full > run.log 2>&1 &

# 调整参数后，快速重绘
python genome_circos_qva_fast.py --mode plot > run.log 2>&1 &

# 自动模式（推荐）
python genome_circos_qva_fast.py > run.log 2>&1 &
"""

import os
import glob
import json
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import colorsys
from pycirclize import Circos

# 配置字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Liberation Sans']
plt.rcParams['axes.unicode_minus'] = False

# 结果输出路径设置
OUTPUT_DIR = "/home/user/xiezheng/genome_circos/genome_circos_test/run"
FIGURE_SIZE = (23, 23)        # 图片大小
FIGURE_DPI = 500              # 分辨率

OUTPUT_IMAGE = f"{OUTPUT_DIR}/circos_comparison.png"
OUTPUT_PDF = f"{OUTPUT_DIR}/circos_comparison.pdf"
OUTPUT_SVG = f"{OUTPUT_DIR}/circos_comparison.svg"
OUTPUT_ANI_CSV = f"{OUTPUT_DIR}/ani_matrix.csv"
CACHE_LINKS = f"{OUTPUT_DIR}/links_cache.json"
CACHE_GENOMES = f"{OUTPUT_DIR}/genomes_info.json"

# ==================== 配置参数 ====================
# 数据配置（修改后需要 mode=full）
REFERENCE_DIR = "reference"
QUERY_DIR = "query"
NUM_QUERY_GENOMES = -1        # 从query目录选择几个基因组（-1表示全部）

# 比对参数（修改后需要 mode=full）
MIN_IDENTITY = 85.0           # 最小相似度阈值（%）
MIN_LENGTH = 500              # 最小区域长度（bp）

# 可视化参数（以下修改后只需 mode=plot，可快速重绘）
SECTOR_SPACE = 0.5            # 基因组块间隔角度
LABEL_RADIUS = 99             # 标签半径
GENOME_TRACK_INNER_QUERY = 92 # query基因组块内径
GENOME_TRACK_INNER_REF = 93   # reference基因组块内径
GENOME_TRACK_OUTER = 98       # 基因组块外径

# 标签样式
LABEL_SIZE_QUERY = 8          # query基因组标签字体大小
LABEL_SIZE_REF = 6            # reference基因组标签字体大小
LABEL_ADJUST_ROTATION = True  # 标签是否径向排列（True=径向，False=竖直）
LABEL_ORIENTATION = "vertical"# 标签方向（vertical/horizontal）

# 连线样式
LINK_ALPHA_BASE_QUERY = 0.4   # query基因组连线透明度
LINK_ALPHA_BASE_REF = 0.25    # reference基因组连线透明度
LINK_WIDTH_QUERY = 0.9        # query基因组连线粗细
LINK_WIDTH_REF = 0.6          # reference基因组连线粗细


def simplify_genome_name(filename):
    """简化基因组名称，提取前两个下划线分隔的部分"""
    basename = os.path.splitext(os.path.basename(filename))[0]
    parts = basename.split('_')
    return '_'.join(parts[:2]) if len(parts) >= 2 else basename


def get_genome_files(directory):
    """
    获取目录中的所有基因组文件
    支持 .fa, .fna, .fasta 扩展名（大小写不敏感）
    """
    files = []
    patterns = ["*.fa", "*.fna", "*.fasta", "*.FA", "*.FNA", "*.FASTA"]
    for pattern in patterns:
        files.extend(glob.glob(os.path.join(directory, pattern)))
    return sorted(set(files))  # 去重并排序


def generate_reference_colors(n):
    """生成参考基因组的冷色调（蓝绿青）"""
    colors = []
    for i in range(n):
        hue = (180 + (90 * i / max(n-1, 1))) / 360
        saturation = 0.6 + (i % 3) * 0.1
        value = 0.75 + (i % 2) * 0.15
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append('#{:02x}{:02x}{:02x}'.format(
            int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)))
    return colors


def generate_query_colors(n):
    """生成查询基因组的暖色调（红橙黄）"""
    colors = []
    for i in range(n):
        hue = (0 + (60 * i / max(n-1, 1))) / 360
        saturation = 0.8 + (i % 2) * 0.1
        value = 0.85 + (i % 2) * 0.1
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append('#{:02x}{:02x}{:02x}'.format(
            int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)))
    return colors


def check_cache_exists():
    """检查缓存文件是否存在"""
    return (os.path.exists(CACHE_LINKS) and 
            os.path.exists(CACHE_GENOMES) and
            os.path.exists(OUTPUT_ANI_CSV))


def save_links_cache(all_links, genomes):
    """保存连线数据和基因组信息到缓存文件"""
    print(f"\nSaving cache files...")
    
    with open(CACHE_LINKS, 'w') as f:
        json.dump(all_links, f, indent=2)
    print(f"  Links cache: {CACHE_LINKS}")
    
    genomes_data = []
    for g in genomes:
        genomes_data.append({
            'name': g['name'],
            'short_name': g['short_name'],
            'file': g['file'],
            'length': g['length'],
            'offsets': g.get('offsets', {}),  # 保存contig偏移量
            'color': g['color'],
            'is_reference': g['is_reference']
        })
    
    with open(CACHE_GENOMES, 'w') as f:
        json.dump(genomes_data, f, indent=2)
    print(f"  Genomes info: {CACHE_GENOMES}")


def load_cache():
    """从缓存加载数据"""
    print(f"\nLoading cache files...")
    
    with open(CACHE_LINKS, 'r') as f:
        all_links = json.load(f)
    print(f"  Loaded {len(all_links)} link groups")
    
    with open(CACHE_GENOMES, 'r') as f:
        genomes = json.load(f)
    print(f"  Loaded {len(genomes)} genomes")
    
    return all_links, genomes


def check_tools():
    """检查必需工具是否可用"""
    tools = {}
    for tool in ['nucmer', 'fastANI']:
        try:
            cmd = [tool, '--version'] if tool == 'nucmer' else [tool, '-h']
            result = subprocess.run(cmd, capture_output=True, timeout=5)
            tools[tool] = result.returncode == 0 or tool == 'fastANI'
        except:
            tools[tool] = False
    return tools


def calculate_ani(query_file, ref_file):
    """计算两个基因组之间的ANI值"""
    try:
        output = "temp_ani.txt"
        cmd = ['fastANI', '-q', query_file, '-r', ref_file, '-o', output]
        subprocess.run(cmd, capture_output=True, timeout=300, check=False)
        
        if os.path.exists(output):
            with open(output) as f:
                line = f.readline().strip()
                if line:
                    ani = float(line.split('\t')[2])
                    os.remove(output)
                    return ani
        return None
    except:
        return None


def calculate_all_ani(files):
    """计算所有基因组对之间的ANI矩阵"""
    n = len(files)
    names = [simplify_genome_name(f) for f in files]
    ani_matrix = pd.DataFrame(np.nan, index=names, columns=names)
    
    # 对角线设为100%
    for i in range(n):
        ani_matrix.loc[names[i], names[i]] = 100.0
    
    print("\n" + "="*60)
    print("Computing ANI values")
    print("="*60)
    
    total = n * (n - 1) // 2
    current = 0
    success_count = 0
    
    for i in range(n):
        for j in range(i+1, n):
            current += 1
            print(f"[{current}/{total}] {names[i]} vs {names[j]}")
            
            ani = calculate_ani(files[i], files[j])
            if ani:
                ani_matrix.loc[names[i], names[j]] = ani
                ani_matrix.loc[names[j], names[i]] = ani
                print(f"         ANI = {ani:.2f}%")
                success_count += 1
            else:
                print(f"         ANI = Failed (NaN)")
    
    print(f"\nSuccessful calculations: {success_count}/{total}")
    return ani_matrix


def cleanup_nucmer_temp(prefix):
    """
    清理nucmer产生的临时文件
    nucmer会产生多个中间文件，运行结束后需要清理
    """
    extensions = ['.delta', '.ntref', '.mgaps', '.seq', '.coords']
    for ext in extensions:
        temp_file = f"{prefix}{ext}"
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
            except:
                pass


def run_nucmer(file1, file2, prefix):
    """运行nucmer比对"""
    try:
        cmd = ['nucmer', '--maxmatch', '-c', '100', '-p', prefix, file1, file2]
        subprocess.run(cmd, capture_output=True, timeout=600, check=False)
        delta = f"{prefix}.delta"
        return delta if os.path.exists(delta) else None
    except:
        return None


def get_genome_info(fasta_file):
    """
    获取基因组的总长度及每个contig的全局偏移量
    这是修复多contig基因组坐标错误的关键函数
    
    返回:
        total_length: 基因组总长度
        offsets: 字典，contig_id -> 该contig在全局坐标系中的起始位置
    """
    total_length = 0
    offsets = {}
    current_seq_id = None
    current_len = 0

    try:
        with open(fasta_file) as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 保存上一个contig的信息
                    if current_seq_id is not None:
                        offsets[current_seq_id] = total_length
                        total_length += current_len
                    
                    # 提取新的contig ID（取>后第一个单词）
                    current_seq_id = line[1:].split()[0]
                    current_len = 0
                else:
                    current_len += len(line)
            
            # 保存最后一个contig
            if current_seq_id is not None:
                offsets[current_seq_id] = total_length
                total_length += current_len
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}")
        return 5000000, {}
    
    # 验证结果
    if total_length == 0:
        print(f"Warning: {fasta_file} appears to be empty")
        return 5000000, {}
    
    if len(offsets) == 0:
        print(f"Warning: No sequences found in {fasta_file}")
        return 5000000, {}

    return total_length, offsets


def parse_delta(delta_file, min_id, min_len, ref_offsets, query_offsets):
    """
    解析nucmer的delta文件，将局部坐标转换为全局坐标
    
    关键修复：
    nucmer输出的是相对于每个contig的局部坐标
    Circos需要的是相对于整个基因组的全局坐标
    必须加上contig的偏移量才能得到正确的位置
    
    参数:
        delta_file: nucmer输出的delta文件路径
        min_id: 最小相似度阈值
        min_len: 最小区域长度
        ref_offsets: 参考基因组的contig偏移量字典
        query_offsets: 查询基因组的contig偏移量字典
    
    返回:
        links: 相似区域列表，每个元素包含全局坐标
    """
    links = []
    current_ref_id = None
    current_query_id = None

    try:
        with open(delta_file) as f:
            for line in f:
                line = line.strip()
                
                # 跳过文件头
                if not line or line.startswith('/') or line in ['NUCMER', 'PROMER']:
                    continue
                
                # 捕获当前正在比对的contig名称
                if line.startswith('>'):
                    parts = line[1:].split()
                    current_ref_id = parts[0]
                    current_query_id = parts[1]
                    continue
                
                # 解析比对坐标
                parts = line.split()
                if len(parts) == 7:
                    r1, r2, q1, q2 = map(int, parts[:4])
                    errors = int(parts[4])
                    length = max(abs(r2 - r1), abs(q2 - q1))
                    
                    if length > 0:
                        # 计算相似度
                        identity = (1 - errors / length) * 100
                        
                        # 过滤：只保留满足阈值的区域
                        if identity >= min_id and length >= min_len:
                            # 获取当前contig的偏移量
                            r_off = ref_offsets.get(current_ref_id, 0)
                            q_off = query_offsets.get(current_query_id, 0)
                            
                            # 关键步骤：局部坐标 + 偏移量 = 全局坐标
                            links.append({
                                'ref_start': min(r1, r2) + r_off,
                                'ref_end': max(r1, r2) + r_off,
                                'query_start': min(q1, q2) + q_off,
                                'query_end': max(q1, q2) + q_off,
                                'identity': identity,
                                'length': length
                            })
    except Exception as e:
        print(f"    Warning parsing delta: {e}")
    
    return links


def compute_alignments(genomes):
    """
    执行Query vs All的全基因组比对策略
    
    策略说明：
    只有query基因组会与所有基因组（包括reference和其他query）进行比对
    这样避免了不必要的reference vs reference比对
    """
    print("\n" + "="*60)
    print("Performing genome alignments (Query vs All)")
    print("="*60)
    
    query_genomes = [g for g in genomes if not g["is_reference"]]
    all_genomes = genomes
    
    # 构建比对配对列表
    pairs = []
    for q in query_genomes:
        for other in all_genomes:
            if q["short_name"] != other["short_name"]:
                pairs.append((q, other))
    
    print(f"Query genomes: {len(query_genomes)}")
    print(f"Total genomes: {len(all_genomes)}")
    print(f"Total comparisons: {len(pairs)} (Query vs All)")
    print("-" * 60)
    
    all_links = []
    
    for idx, (g1, g2) in enumerate(pairs, 1):
        prefix = f"temp_{idx}"
        
        g1_type = "Query" if not g1["is_reference"] else "Ref"
        g2_type = "Query" if not g2["is_reference"] else "Ref"
        print(f"[{idx}/{len(pairs)}] [{g1_type}] {g1['short_name']} <-> [{g2_type}] {g2['short_name']}")
        
        # 运行nucmer比对（g1作为参考，g2作为查询）
        delta = run_nucmer(g1["file"], g2["file"], prefix)
        
        if delta:
            # 解析delta文件，应用坐标转换
            links = parse_delta(delta, MIN_IDENTITY, MIN_LENGTH, 
                              g1['offsets'], g2['offsets'])
            print(f"  -> Found {len(links)} similar regions")
            
            if links:
                all_links.append({
                    'g1_name': g1['short_name'],
                    'g2_name': g2['short_name'],
                    'g1_is_ref': g1['is_reference'],
                    'links': links
                })
        else:
            print(f"  -> Alignment failed")
        
        # 清理本轮比对产生的临时文件
        cleanup_nucmer_temp(prefix)
    
    print("-" * 60)
    print(f"Total link groups: {len(all_links)}")
    
    return all_links


def plot_circos(genomes, all_links):
    """绘制Circos图"""
    print("\n" + "="*60)
    print("Plotting Circos diagram")
    print("="*60)
    
    # 构建sectors（每个基因组一个扇区）
    sectors = {g['short_name']: g['length'] for g in genomes}
    
    for g in genomes:
        genome_type = "Query" if not g['is_reference'] else "Ref"
        print(f"  [{genome_type}] {g['short_name']}: {g['length']:,} bp")
    
    print(f"\nVisualization parameters:")
    print(f"  Sector space: {SECTOR_SPACE} degrees")
    print(f"  Label radius: {LABEL_RADIUS}")
    print(f"  Label rotation: {'Radial' if LABEL_ADJUST_ROTATION else 'Fixed'}")
    print(f"  Label orientation: {LABEL_ORIENTATION if not LABEL_ADJUST_ROTATION else 'Auto'}")
    print(f"  Figure size: {FIGURE_SIZE}")
    print(f"  DPI: {FIGURE_DPI}")
    
    circos = Circos(sectors, space=SECTOR_SPACE)
    
    # 绘制基因组轨道
    print("\nDrawing genome tracks...")
    for sector in circos.sectors:
        g = next(g for g in genomes if g["short_name"] == sector.name)
        
        # query和reference使用不同的轨道内径
        track_inner = GENOME_TRACK_INNER_QUERY if not g["is_reference"] else GENOME_TRACK_INNER_REF
        track = sector.add_track((track_inner, GENOME_TRACK_OUTER))
        
        # query和reference使用不同的线宽
        linewidth = 2.5 if not g["is_reference"] else 1.5
        track.axis(fc=g["color"], lw=linewidth, ec='black')
        
        # 添加基因组名称标签
        label_size = LABEL_SIZE_QUERY if not g["is_reference"] else LABEL_SIZE_REF
        label_weight = 'bold' if not g["is_reference"] else 'normal'
        
        # 构建标签参数
        sector.text(
            text=sector.name,
            r=LABEL_RADIUS,
            size=label_size,
            orientation=LABEL_ORIENTATION,         # 使用开头的LABEL_ORIENTATION设置
            adjust_rotation=LABEL_ADJUST_ROTATION, # 使用开头的LABEL_ADJUST_ROTATION设置
            fontweight=label_weight
        )
    
    # 绘制连线
    print("Drawing similarity links...")
    total_links_drawn = 0
    genome_map = {g['short_name']: g for g in genomes}
    
    for link_group in all_links:
        g1 = genome_map[link_group['g1_name']]
        g2 = genome_map[link_group['g2_name']]
        
        for link in link_group['links']:
            try:
                # 根据相似度计算透明度
                alpha_base = LINK_ALPHA_BASE_QUERY if not g1["is_reference"] else LINK_ALPHA_BASE_REF
                alpha_bonus = (link['identity'] - MIN_IDENTITY) / 50 * 0.2
                alpha = min(0.6, alpha_base + alpha_bonus)
                
                # 设置线宽
                linewidth = LINK_WIDTH_QUERY if not g1["is_reference"] else LINK_WIDTH_REF
                
                # 绘制连线
                # 注意：nucmer中g1是ref，g2是query
                # 所以g1对应ref_start/end，g2对应query_start/end
                circos.link(
                    (g1["short_name"], link['ref_start'], link['ref_end']),
                    (g2["short_name"], link['query_start'], link['query_end']),
                    color=g1["color"],
                    alpha=alpha,
                    lw=linewidth
                )
                total_links_drawn += 1
            except:
                pass
    
    print(f"Total links drawn: {total_links_drawn}")
    
    # 渲染图形
    print(f"\nRendering figure ({FIGURE_SIZE[0]}x{FIGURE_SIZE[1]} @ {FIGURE_DPI} DPI)...")
    fig = circos.plotfig(figsize=FIGURE_SIZE)
    
    # 添加图例
    ref_count = sum(1 for g in genomes if g['is_reference'])
    query_count = len(genomes) - ref_count
    
    ref_patch = mpatches.Patch(color='#4A90E2', 
                               label=f'Reference Genomes (n={ref_count})')
    query_patch = mpatches.Patch(color='#E74C3C', 
                                 label=f'Query Genomes (n={query_count})')
    
    fig.legend(
        handles=[ref_patch, query_patch],
        loc='lower center',
        ncol=2,
        frameon=True,
        fancybox=True,
        shadow=True,
        fontsize=14,
        title='Genome Types',
        title_fontsize=16,
        bbox_to_anchor=(0.5, -0.02)
    )
    
    plt.subplots_adjust(left=0.15, right=0.85, top=0.90, bottom=0.05)
    
    # 保存图片
    print(f"\nSaving images...")
    fig.savefig(OUTPUT_IMAGE, dpi=FIGURE_DPI, bbox_inches='tight')
    print(f"  PNG: {OUTPUT_IMAGE}")
    
    fig.savefig(OUTPUT_PDF, bbox_inches='tight')
    print(f"  PDF: {OUTPUT_PDF}")
    
    fig.savefig(OUTPUT_SVG, format='svg', bbox_inches='tight')
    print(f"  SVG: {OUTPUT_SVG}")
    
    plt.close()
    print("\nCircos plot completed!")


def main():
    parser = argparse.ArgumentParser(
        description='Genome Circos Visualization - Query vs All with Fast Replot')
    parser.add_argument('--mode', 
                       choices=['full', 'plot', 'auto'], 
                       default='auto',
                       help='full: complete run | plot: plot only | auto: auto detect (default)')
    args = parser.parse_args()
    
    # 创建输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print("=" * 60)
    print("Genome Circos - Query vs All - Fast Replot Edition")
    print("=" * 60)
    
    # 检查缓存
    cache_exists = check_cache_exists()
    
    if args.mode == 'auto':
        mode = 'plot' if cache_exists else 'full'
        print(f"\nAuto mode: Detected {'cache files' if cache_exists else 'no cache'}")
        print(f"  Running in '{mode}' mode")
    else:
        mode = args.mode
        if mode == 'plot' and not cache_exists:
            print("\nError: plot mode requires cache files")
            print("  Please run with --mode full first")
            return
    
    print(f"\nRun mode: {mode.upper()}")
    print(f"Comparison strategy: Query vs All")
    
    if mode == 'plot':
        # 快速重绘模式：只读取缓存并绘图
        print("\nFast replot mode - only visualization will be updated")
        print("  Estimated time: ~30 seconds")
        
        all_links, genomes = load_cache()
        plot_circos(genomes, all_links)
        
    else:
        # 完整运行模式
        print("\nFull computation mode")
        print("  Future runs can use --mode plot for instant results")
        
        # 检查必需工具
        print("\nChecking tools...")
        tools = check_tools()
        print(f"  nucmer: {'Available' if tools['nucmer'] else 'Not found'}")
        print(f"  fastANI: {'Available' if tools['fastANI'] else 'Not found'}")
        
        if not tools['nucmer']:
            print("\nError: nucmer is required but not found")
            return
        
        # 获取基因组文件列表
        ref_files = get_genome_files(REFERENCE_DIR)
        query_files = get_genome_files(QUERY_DIR)
        
        if not ref_files:
            print(f"\nError: No genome files found in {REFERENCE_DIR}")
            return
        
        if not query_files:
            print(f"\nError: No genome files found in {QUERY_DIR}")
            return
        
        print(f"\nFound {len(ref_files)} reference genomes")
        print(f"Found {len(query_files)} query genomes")
        
        # 限制query基因组数量（如果设置了）
        if NUM_QUERY_GENOMES > 0:
            query_files = query_files[:NUM_QUERY_GENOMES]
            print(f"Selected {len(query_files)} query genomes")
        
        # 计算ANI矩阵（如果fastANI可用）
        all_files = ref_files + query_files
        if tools['fastANI']:
            ani_matrix = calculate_all_ani(all_files)
            ani_matrix.to_csv(OUTPUT_ANI_CSV)
            print(f"\nANI matrix saved: {OUTPUT_ANI_CSV}")
        else:
            print("\nSkipping ANI calculation (fastANI not available)")
        
        # 生成颜色方案
        ref_colors = generate_reference_colors(len(ref_files))
        query_colors = generate_query_colors(len(query_files))
        
        # 处理参考基因组
        genomes = []
        for i, f in enumerate(ref_files):
            short_name = simplify_genome_name(f)
            length, offsets = get_genome_info(f)
            genomes.append({
                "name": short_name,
                "short_name": short_name,
                "file": f,
                "length": length,
                "offsets": offsets,
                "color": ref_colors[i],
                "is_reference": True
            })
        
        # 处理查询基因组
        for i, f in enumerate(query_files):
            short_name = simplify_genome_name(f)
            length, offsets = get_genome_info(f)
            genomes.append({
                "name": short_name,
                "short_name": short_name,
                "file": f,
                "length": length,
                "offsets": offsets,
                "color": query_colors[i],
                "is_reference": False
            })
        
        # 执行比对
        all_links = compute_alignments(genomes)
        
        # 保存缓存
        save_links_cache(all_links, genomes)
        
        # 绘制图形
        plot_circos(genomes, all_links)
    
    print("\n" + "=" * 60)
    print("All tasks completed!")
    print("=" * 60)
    print(f"\nTip: To adjust visual parameters and replot instantly:")
    print(f"  1. Edit visualization parameters in the script")
    print(f"  2. Run: python {os.path.basename(__file__)} --mode plot")
    print(f"  3. Done in ~30 seconds!")
    print()


if __name__ == "__main__":
    main()