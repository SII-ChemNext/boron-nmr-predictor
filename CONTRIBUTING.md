# 贡献指南

感谢您对硼核磁预测系统项目的关注！我们欢迎各种形式的贡献。

## 如何贡献

### 报告 Bug

如果您发现了 bug，请：

1. 检查 [Issues](https://github.com/your-username/boron-nmr-predictor/issues) 确认问题未被报告
2. 创建新 Issue，包含：
   - 清晰的标题和描述
   - 复现步骤
   - 预期行为和实际行为
   - 系统环境（OS、Python 版本等）
   - 相关截图或错误日志

### 提出新功能

1. 先创建 Issue 讨论新功能的必要性和设计
2. 等待维护者反馈
3. 获得批准后再开始开发

### 提交代码

1. **Fork 项目**
   ```bash
   # 在 GitHub 上点击 Fork 按钮
   git clone https://github.com/your-username/boron-nmr-predictor.git
   cd boron-nmr-predictor
   ```

2. **创建分支**
   ```bash
   git checkout -b feature/your-feature-name
   # 或
   git checkout -b fix/your-bug-fix
   ```

3. **编写代码**
   - 遵循项目代码风格
   - 添加必要的注释
   - 更新相关文档

4. **测试**
   ```bash
   # 运行测试
   pytest tests/
   
   # 检查代码风格
   flake8 web_app/
   black --check web_app/
   ```

5. **提交更改**
   ```bash
   git add .
   git commit -m "feat: 添加新功能描述"
   # 或
   git commit -m "fix: 修复 bug 描述"
   ```

6. **推送到 GitHub**
   ```bash
   git push origin feature/your-feature-name
   ```

7. **创建 Pull Request**
   - 在 GitHub 上创建 PR
   - 填写 PR 模板
   - 等待代码审查

## 代码规范

### Python 代码风格

- 遵循 [PEP 8](https://pep8.org/)
- 使用 Black 格式化代码
- 使用类型提示（Type Hints）
- 函数和类添加文档字符串

示例：
```python
def predict_chemical_shift(
    smiles: str, 
    solvent: str
) -> Dict[str, Any]:
    """
    预测分子的化学位移
    
    Args:
        smiles: 分子的 SMILES 字符串
        solvent: 溶剂名称
        
    Returns:
        包含预测结果的字典
        
    Raises:
        ValueError: 当 SMILES 无效时
    """
    pass
```

### JavaScript 代码风格

- 使用 ES6+ 语法
- 使用 const/let 而非 var
- 添加适当的注释
- 函数名使用驼峰命名

### HTML/CSS 规范

- 语义化 HTML
- 响应式设计
- 可访问性（ARIA 标签）
- CSS 类名使用 kebab-case

## 提交信息规范

使用语义化提交信息：

- `feat:` 新功能
- `fix:` Bug 修复
- `docs:` 文档更新
- `style:` 代码格式（不影响功能）
- `refactor:` 重构
- `test:` 测试相关
- `chore:` 构建/工具相关

示例：
```
feat: 添加批量预测功能
fix: 修复溶剂选择器的显示问题
docs: 更新 API 文档
```

## 开发环境设置

```bash
# 克隆项目
git clone https://github.com/your-username/boron-nmr-predictor.git
cd boron-nmr-predictor/web_app

# 创建虚拟环境
python -m venv venv
source venv/bin/activate

# 安装依赖（包括开发依赖）
pip install -r requirements.txt
pip install pytest pytest-cov black flake8 mypy

# 运行测试
pytest tests/

# 启动开发服务器
python app.py
```

## 测试

### 编写测试

- 为新功能添加单元测试
- 确保测试覆盖率 > 80%
- 测试文件放在 `tests/` 目录

示例：
```python
# tests/test_predictor.py
def test_predict_valid_smiles():
    predictor = BoronNMRPredictor()
    result = predictor.predict("OB(O)c1ccccc1", "CDCl3")
    assert result is not None
    assert len(result) > 0
```

### 运行测试

```bash
# 运行所有测试
pytest

# 运行特定测试
pytest tests/test_predictor.py

# 生成覆盖率报告
pytest --cov=web_app tests/
```

## 文档

- 更新 README.md（如果添加新功能）
- 更新 API 文档
- 添加代码注释
- 更新 CHANGELOG.md

## 代码审查

所有 PR 都需要经过代码审查：

- 至少一位维护者批准
- 所有测试通过
- 代码风格检查通过
- 无合并冲突

## 行为准则

- 尊重所有贡献者
- 建设性的反馈
- 专注于技术讨论
- 欢迎新手

## 许可证

提交代码即表示您同意将代码以 MIT 许可证发布。

## 问题？

如有疑问，请：
- 查看现有 Issues
- 创建新 Issue
- 联系维护者

感谢您的贡献！🎉
