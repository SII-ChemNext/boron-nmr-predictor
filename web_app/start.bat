@echo off
REM 硼核磁预测系统 Web UI 启动脚本 (Windows)

setlocal enabledelayedexpansion

echo ==================================================
echo 硼核磁预测系统 Web UI - 启动脚本
echo ==================================================
echo.

REM 获取当前目录
set SCRIPT_DIR=%cd%
echo 工作目录: %SCRIPT_DIR%

REM 检查 Python
echo 检查 Python 版本...
python --version >nul 2>&1
if errorlevel 1 (
    echo ✗ 错误: 未找到 Python
    pause
    exit /b 1
)

for /f "tokens=*" %%i in ('python --version 2^>^&1') do set PYTHON_VERSION=%%i
echo ✓ %PYTHON_VERSION%
echo.

REM 检查虚拟环境
if not exist "venv" (
    echo 创建虚拟环境...
    python -m venv venv
    echo ✓ 虚拟环境已创建
)

REM 激活虚拟环境
echo 激活虚拟环境...
call venv\Scripts\activate.bat
echo ✓ 虚拟环境已激活
echo.

REM 检查依赖
echo 检查依赖包...
if not exist "requirements.txt" (
    echo ✗ 错误: 找不到 requirements.txt
    pause
    exit /b 1
)

pip install -q -r requirements.txt
echo ✓ 依赖包已安装
echo.

REM 检查模型文件
echo 检查模型文件...
set MISSING_MODELS=0
for /L %%i in (1,1,5) do (
    if not exist "models\model_fold_%%i.pth" (
        echo ✗ 缺少模型文件: models\model_fold_%%i.pth
        set MISSING_MODELS=1
    )
)

if exist "core\model.py" if exist "core\features.py" (
    echo ✓ 模型文件检查完成
) else (
    echo ✗ 错误: 缺少必要的模型文件
    pause
    exit /b 1
)
echo.

REM 创建必要的目录
echo 创建必要的目录...
if not exist "database" mkdir database
if not exist "static\img" mkdir static\img
echo ✓ 目录已创建
echo.

REM 检查数据库
if not exist "database\predictions.db" (
    echo 初始化数据库...
    python -c "from database.models import init_db; init_db('database\predictions.db')"
    echo ✓ 数据库已初始化
) else (
    echo ✓ 数据库已存在
)
echo.

REM 启动应用
echo ==================================================
echo 启动应用...
echo ==================================================
echo.
echo 服务地址: http://localhost:5000
echo.
echo 按 Ctrl+C 停止应用
echo.

if "%1"=="production" (
    echo 使用生产模式 (Gunicorn) 启动...
    pip install -q gunicorn
    gunicorn -w 4 -b 0.0.0.0:5000 --timeout 120 app:app
) else (
    echo 使用开发模式 (Flask) 启动...
    python app.py
)

pause
