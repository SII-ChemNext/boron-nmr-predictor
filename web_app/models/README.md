# 模型文件

## 下载说明

由于模型文件较大（约 500MB），未包含在 Git 仓库中。

请从以下链接下载模型文件：

- **Google Drive**: https://drive.google.com/drive/folders/1HLirbH9JOf6HvgwUTIUETkpXG_U4fPYl?usp=sharing

## 文件列表

下载后，请将以下 5 个模型文件放置在此目录：

```
models/
├── model_fold_1.pth  (~100MB)
├── model_fold_2.pth  (~100MB)
├── model_fold_3.pth  (~100MB)
├── model_fold_4.pth  (~100MB)
└── model_fold_5.pth  (~100MB)
```

## 模型信息

- **架构**: Graph Neural Network (GNN)
- **训练方法**: 5-Fold Cross Validation
- **输入**: 分子图 + 溶剂指纹
- **输出**: 硼原子化学位移 (ppm)
- **参数量**: 约 2M per model

## 验证

模型文件下载完成后，可以运行以下命令验证：

```bash
python -c "import torch; print('Model 1:', torch.load('model_fold_1.pth', map_location='cpu').keys())"
```

如果输出模型的键值，说明文件正确。
