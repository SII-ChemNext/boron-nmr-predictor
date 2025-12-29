#!/bin/bash

# 硼核磁预测系统 Web UI 启动脚本

set -e

echo "=================================================="
echo "硼核磁预测系统 Web UI - 启动脚本"
echo "=================================================="
echo ""

# 获取当前脚本目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "工作目录: $SCRIPT_DIR"
cd "$SCRIPT_DIR"

# 检查 Python 版本
echo "检查 Python 版本..."
if ! command -v python3 &> /dev/null; then
    echo "✗ 错误: 未找到 Python 3"
    exit 1
fi

PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "✓ Python 版本: $PYTHON_VERSION"
echo ""

# 检查虚拟环境
if [ ! -d "venv" ]; then
    echo "创建虚拟环境..."
    python3 -m venv venv
    echo "✓ 虚拟环境已创建"
fi

# 激活虚拟环境
echo "激活虚拟环境..."
source venv/bin/activate
echo "✓ 虚拟环境已激活"
echo ""

# 检查依赖
echo "检查依赖包..."
if [ ! -f "requirements.txt" ]; then
    echo "✗ 错误: 找不到 requirements.txt"
    exit 1
fi

pip install -q -r requirements.txt
echo "✓ 依赖包已安装"
echo ""

# 检查模型文件
echo "检查模型文件..."
MISSING_MODELS=0
for i in {1..5}; do
    if [ ! -f "models/model_fold_$i.pth" ]; then
        echo "✗ 缺少模型文件: models/model_fold_$i.pth"
        MISSING_MODELS=1
    fi
done

if [ -f "core/model.py" ] && [ -f "core/features.py" ]; then
    echo "✓ 模型文件检查完成"
else
    echo "✗ 错误: 缺少必要的模型文件"
    exit 1
fi
echo ""

# 创建必要的目录
echo "创建必要的目录..."
mkdir -p database static/img
echo "✓ 目录已创建"
echo ""

# 检查数据库
if [ ! -f "database/predictions.db" ]; then
    echo "初始化数据库..."
    python3 << EOF
from database.models import init_db
init_db('database/predictions.db')
EOF
    echo "✓ 数据库已初始化"
else
    echo "✓ 数据库已存在"
fi
echo ""

# 启动应用
echo "=================================================="
echo "启动应用..."
echo "=================================================="
echo ""
echo "服务地址: http://localhost:5000"
echo ""
echo "按 Ctrl+C 停止应用"
echo ""

# 检查是否使用 Gunicorn（生产环境）
if [ "$1" == "production" ]; then
    echo "使用生产模式 (Gunicorn) 启动..."
    pip install -q gunicorn
    gunicorn -w 4 -b 0.0.0.0:5000 --timeout 120 app:app
else
    echo "使用开发模式 (Flask) 启动..."
    python3 app.py
fi
