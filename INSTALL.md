# 安装指南

## 快速安装（5 分钟）

### 1. 克隆项目

```bash
git clone https://github.com/your-username/boron-nmr-predictor.git
cd boron-nmr-predictor/web_app
```

### 2. 创建虚拟环境

**Linux / macOS:**
```bash
python3 -m venv venv
source venv/bin/activate
```

**Windows:**
```bash
python -m venv venv
venv\Scripts\activate
```

### 3. 安装依赖

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

**注意**：安装 PyTorch 和 RDKit 可能需要 5-10 分钟，请耐心等待。

### 4. 下载模型文件

从以下链接下载模型文件（约 500MB）：

- **Google Drive**: https://drive.google.com/drive/folders/1HLirbH9JOf6HvgwUTIUETkpXG_U4fPYl?usp=sharing

将下载的 5 个 `.pth` 文件放入 `models/` 目录：

```bash
# 确保文件在正确位置
ls models/
# 应该看到：
# model_fold_1.pth  model_fold_2.pth  model_fold_3.pth  
# model_fold_4.pth  model_fold_5.pth
```

### 5. 启动应用

**Linux / macOS:**
```bash
chmod +x start.sh
./start.sh
```

**Windows:**
```bash
start.bat
```

**或手动启动:**
```bash
python app.py
```

### 6. 访问应用

在浏览器中打开：http://localhost:5000

## 详细安装说明

### 系统要求

- **操作系统**: Linux / macOS / Windows
- **Python**: 3.8 或更高版本
- **内存**: 最低 4GB，推荐 8GB
- **磁盘空间**: 约 2GB（包含依赖和模型）
- **网络**: 首次安装需要下载依赖包

### 可选：GPU 加速

如果您有 NVIDIA GPU，可以安装 CUDA 版本的 PyTorch 以加速预测：

```bash
# 卸载 CPU 版本
pip uninstall torch

# 安装 CUDA 版本（以 CUDA 11.8 为例）
pip install torch==2.1.0 --index-url https://download.pytorch.org/whl/cu118
```

### 依赖说明

主要依赖包：

- **Flask 3.0.0** - Web 框架
- **PyTorch 2.1.0** - 深度学习框架
- **PyTorch Geometric 2.4.0** - 图神经网络
- **RDKit 2023.9.1** - 化学信息学
- **NumPy, Pandas** - 数据处理
- **Pillow, Matplotlib** - 图像处理

完整依赖列表见 `requirements.txt`

### 常见安装问题

#### 问题 1: RDKit 安装失败

**解决方案 1 - 使用 conda:**
```bash
conda install -c conda-forge rdkit
```

**解决方案 2 - 使用预编译包:**
```bash
pip install rdkit-pypi
```

#### 问题 2: PyTorch Geometric 安装失败

**解决方案:**
```bash
pip install torch-scatter torch-sparse torch-cluster torch-spline-conv -f https://data.pyg.org/whl/torch-2.1.0+cpu.html
pip install torch-geometric
```

#### 问题 3: 内存不足

如果在安装过程中遇到内存不足：

```bash
# 逐个安装大型包
pip install torch
pip install torch-geometric
pip install rdkit
pip install -r requirements.txt
```

#### 问题 4: 权限错误（Linux/Mac）

```bash
# 使用 --user 标志
pip install --user -r requirements.txt
```

### 验证安装

运行以下命令验证安装是否成功：

```bash
python -c "import torch; import torch_geometric; import rdkit; print('✓ 所有依赖安装成功')"
```

### 开发环境设置

如果您想参与开发，还需要安装开发依赖：

```bash
pip install pytest pytest-cov black flake8 mypy
```

## 卸载

```bash
# 停用虚拟环境
deactivate

# 删除项目目录
cd ..
rm -rf boron-nmr-predictor
```

## 更新

```bash
cd boron-nmr-predictor
git pull origin main
pip install -r web_app/requirements.txt --upgrade
```

## 获取帮助

如果遇到问题：

1. 查看 [常见问题](README.md#常见问题)
2. 搜索 [Issues](https://github.com/your-username/boron-nmr-predictor/issues)
3. 提交新的 Issue

## 下一步

安装完成后，请阅读 [README.md](README.md) 了解如何使用应用。
