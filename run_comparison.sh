#!/bin/bash
# 橙皮苷逆合成分析比较脚本
# 分别运行baseline、flavonoid和优化配置，然后比较结果

set -e  # 遇到错误立即退出

echo "=========================================="
echo "橙皮苷逆合成分析比较"
echo "=========================================="
echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# 激活conda环境
source /home/ljx/miniforge3/etc/profile.d/conda.sh
conda activate retro

# 切换到项目目录（flavonoid配置需要加载scripts模块）
cd /home/ljx/retro_synthesis
export PYTHONPATH="/home/ljx/retro_synthesis:$PYTHONPATH"

# 设置路径
SCRIPTS_DIR="/home/ljx/retro_synthesis/scripts"
OUTPUT_DIR="/home/ljx/retro_synthesis/outputs"

# 确保输出目录存在
mkdir -p "$OUTPUT_DIR"

# 运行baseline分析
echo "=========================================="
echo "第1步: 运行Baseline配置（仅USPTO模板）"
echo "=========================================="
python "$SCRIPTS_DIR/run_baseline_hesperidin.py"

echo ""
echo "=========================================="
echo "第2步: 运行Flavonoid配置（USPTO + 黄酮类模板）"
echo "=========================================="
python "$SCRIPTS_DIR/run_flavonoid_hesperidin.py"

echo ""
echo "=========================================="
echo "第3步: 运行优化配置（改进的MCTS参数）"
echo "=========================================="
python "$SCRIPTS_DIR/run_optimized_hesperidin.py"

# 查找最新的输出文件
echo ""
echo "=========================================="
echo "第4步: 比较结果"
echo "=========================================="

# 获取今天生成的文件
TODAY=$(date '+%y-%m-%d')
BASELINE_FILE=$(ls -t "$OUTPUT_DIR"/${TODAY}-hesperidin_baseline_uspto.json 2>/dev/null | head -1)
FLAVONOID_FILE=$(ls -t "$OUTPUT_DIR"/${TODAY}-hesperidin_flavonoid_combined.json 2>/dev/null | head -1)
OPTIMIZED_FILE=$(ls -t "$OUTPUT_DIR"/${TODAY}-hesperidin_optimized.json 2>/dev/null | head -1)

if [ -z "$BASELINE_FILE" ] || [ -z "$FLAVONOID_FILE" ] || [ -z "$OPTIMIZED_FILE" ]; then
    echo "错误: 找不到输出文件"
    echo "Baseline文件: $BASELINE_FILE"
    echo "Flavonoid文件: $FLAVONOID_FILE"
    echo "Optimized文件: $OPTIMIZED_FILE"
    exit 1
fi

echo "Baseline文件: $BASELINE_FILE"
echo "Flavonoid文件: $FLAVONOID_FILE"
echo "Optimized文件: $OPTIMIZED_FILE"

# 运行比较脚本
python "$SCRIPTS_DIR/compare_all_results.py" "$BASELINE_FILE" "$FLAVONOID_FILE" "$OPTIMIZED_FILE"

echo ""
echo "=========================================="
echo "分析完成!"
echo "=========================================="
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo ""
echo "输出文件位置:"
echo "  - Baseline结果: $BASELINE_FILE"
echo "  - Flavonoid结果: $FLAVONOID_FILE"
echo "  - Optimized结果: $OPTIMIZED_FILE"
