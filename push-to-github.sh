#!/bin/bash

# GitHub 推送脚本
# 使用方法：./push-to-github.sh your-username

if [ -z "$1" ]; then
    echo "❌ 错误：请提供 GitHub 用户名"
    echo "使用方法：./push-to-github.sh your-username"
    exit 1
fi

USERNAME=$1
REPO_NAME="boron-nmr-predictor"

echo "╔════════════════════════════════════════════════════════════╗"
echo "║                                                            ║"
echo "║          🚀 推送到 GitHub                                   ║"
echo "║                                                            ║"
echo "╚════════════════════════════════════════════════════════════╝"
echo ""
echo "📦 仓库：$USERNAME/$REPO_NAME"
echo ""

# 检查是否已有 remote
if git remote | grep -q "origin"; then
    echo "⚠️  检测到已存在的 remote，删除旧的..."
    git remote remove origin
fi

# 添加 remote
echo "📡 添加远程仓库..."
git remote add origin https://github.com/$USERNAME/$REPO_NAME.git

# 重命名分支
echo "🔄 重命名分支为 main..."
git branch -M main

# 推送
echo "⬆️  推送到 GitHub..."
echo ""
echo "⚠️  提示：当要求输入密码时，请使用 Personal Access Token"
echo "   不是你的 GitHub 密码！"
echo ""
echo "如何获取 Token："
echo "1. 访问：https://github.com/settings/tokens"
echo "2. 点击 'Generate new token (classic)'"
echo "3. 勾选 'repo' 权限"
echo "4. 复制生成的 Token"
echo ""
read -p "按 Enter 继续推送..."

git push -u origin main

if [ $? -eq 0 ]; then
    echo ""
    echo "╔════════════════════════════════════════════════════════════╗"
    echo "║                                                            ║"
    echo "║          ✅ 推送成功！                                      ║"
    echo "║                                                            ║"
    echo "╚════════════════════════════════════════════════════════════╝"
    echo ""
    echo "🌐 访问你的项目："
    echo "   https://github.com/$USERNAME/$REPO_NAME"
    echo ""
    echo "📝 下一步："
    echo "   1. 添加截图到 docs/images/web-interface.png"
    echo "   2. 在 GitHub 添加 Topics 标签"
    echo "   3. 完善 About 部分"
else
    echo ""
    echo "❌ 推送失败！"
    echo ""
    echo "常见问题："
    echo "1. 确保已在 GitHub 创建仓库"
    echo "2. 使用 Personal Access Token 而不是密码"
    echo "3. 检查网络连接"
fi
