/**
 * 预测结果处理和展示
 */

let currentPredictionData = null;
let currentPredictionId = null;

/**
 * 显示加载状态
 */
function showLoadingState() {
    document.getElementById('empty-state').classList.add('hidden');
    document.getElementById('result-state').classList.add('hidden');
    document.getElementById('error-state').classList.add('hidden');
    document.getElementById('loading-state').classList.remove('hidden');
}

/**
 * 显示结果状态
 */
function showResultState() {
    document.getElementById('empty-state').classList.add('hidden');
    document.getElementById('loading-state').classList.add('hidden');
    document.getElementById('error-state').classList.add('hidden');
    document.getElementById('result-state').classList.remove('hidden');
}

/**
 * 显示错误状态
 *
 * @param {string} message - 错误消息
 */
function showErrorState(message) {
    document.getElementById('empty-state').classList.add('hidden');
    document.getElementById('loading-state').classList.add('hidden');
    document.getElementById('result-state').classList.add('hidden');
    document.getElementById('error-state').classList.remove('hidden');
    document.getElementById('error-message').textContent = message;
}

/**
 * 显示空状态
 */
function showEmptyState() {
    document.getElementById('loading-state').classList.add('hidden');
    document.getElementById('result-state').classList.add('hidden');
    document.getElementById('error-state').classList.add('hidden');
    document.getElementById('empty-state').classList.remove('hidden');
}

/**
 * 处理预测响应
 *
 * @param {Object} data - 来自 API 的响应数据
 */
function displayPredictionResult(data) {
    if (!data.success) {
        const errorMsg = data.error || (window.i18n ? window.i18n.t('messages.predictionFailed') : '预测失败');
        showErrorState(errorMsg);
        return;
    }

    // 保存当前预测数据
    currentPredictionData = data;
    currentPredictionId = data.prediction_id;

    // 显示标准化 SMILES
    const canonicalSmilesEl = document.getElementById('canonical-smiles');
    canonicalSmilesEl.textContent = data.canonical_smiles;

    // 填充预测表格
    const tbody = document.querySelector('#prediction-tbody');
    tbody.innerHTML = '';

    data.predictions.forEach(pred => {
        const row = document.createElement('tr');
        row.innerHTML = `
            <td>${pred.atom_index}</td>
            <td>${pred.element}</td>
            <td>${pred.ppm.toFixed(2)}</td>
        `;
        tbody.appendChild(row);
    });

    // 显示分子图
    if (data.image_url) {
        const imgEl = document.getElementById('molecule-img');
        imgEl.src = data.image_url;
        imgEl.alt = `分子结构 - ${data.canonical_smiles}`;
    }

    // 显示结果状态
    showResultState();
}

/**
 * 页面加载时的初始化
 */
document.addEventListener('DOMContentLoaded', () => {
    const predictBtn = document.getElementById('predict-btn');
    const closeErrorBtn = document.getElementById('close-error-btn');
    const copySmilesBtn = document.getElementById('copy-smiles-btn');
    const downloadCsvBtn = document.getElementById('download-csv-btn');
    const downloadJsonBtn = document.getElementById('download-json-btn');
    const downloadImageBtn = document.getElementById('download-image-btn');

    // 预测按钮点击事件
    if (predictBtn) {
        predictBtn.addEventListener('click', handlePredict);
    }

    // 关闭错误按钮
    if (closeErrorBtn) {
        closeErrorBtn.addEventListener('click', showEmptyState);
    }

    // 复制 SMILES 按钮
    if (copySmilesBtn) {
        copySmilesBtn.addEventListener('click', () => {
            const smilesText = document.getElementById('canonical-smiles').textContent;
            copyToClipboard(smilesText);
        });
    }

    // 下载按钮
    if (downloadCsvBtn) {
        downloadCsvBtn.addEventListener('click', () => {
            if (currentPredictionId) {
                window.location.href = `/api/download/csv/${currentPredictionId}`;
            }
        });
    }

    if (downloadJsonBtn) {
        downloadJsonBtn.addEventListener('click', () => {
            if (currentPredictionId) {
                window.location.href = `/api/download/json/${currentPredictionId}`;
            }
        });
    }

    if (downloadImageBtn) {
        downloadImageBtn.addEventListener('click', downloadMoleculeImage);
    }

    // 初始显示空状态
    showEmptyState();
});

/**
 * 处理预测请求
 */
async function handlePredict() {
    try {
        // 显示加载状态
        showLoadingState();

        // 获取 SMILES
        const smiles = await getSmilesFromKetcher();

        // 获取选中的溶剂
        const solvent = document.getElementById('solvent-select').value;

        // 发送预测请求
        const response = await fetch('/api/predict', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                molecule_smiles: smiles,
                solvent: solvent
            })
        });

        const data = await response.json();

        // 处理响应
        if (response.ok) {
            displayPredictionResult(data);
        } else {
            const errorMsg = data.error || (window.i18n ? window.i18n.t('messages.predictionFailed') : '预测失败');
            showErrorState(errorMsg);
        }

    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.predictionError') : '预测错误';
        console.error(errorPrefix, error);
        const failedMsg = window.i18n ? window.i18n.t('messages.predictionFailed') : '预测失败';
        showErrorState(`${failedMsg}: ${error.message}`);
    }
}

/**
 * 复制文本到剪贴板
 *
 * @param {string} text - 要复制的文本
 */
async function copyToClipboard(text) {
    try {
        if (navigator.clipboard && window.isSecureContext) {
            // 使用现代 Clipboard API
            await navigator.clipboard.writeText(text);
            const successMsg = window.i18n ? window.i18n.t('messages.copiedToClipboard') : '✓ 已复制到剪贴板';
            alert(successMsg);
        } else {
            // 备选方案：使用传统方法
            const textarea = document.createElement('textarea');
            textarea.value = text;
            document.body.appendChild(textarea);
            textarea.select();
            document.execCommand('copy');
            document.body.removeChild(textarea);
            const successMsg = window.i18n ? window.i18n.t('messages.copiedToClipboard') : '✓ 已复制到剪贴板';
            alert(successMsg);
        }
    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.copyFailed') : '复制失败';
        console.error(errorPrefix, error);
        alert(`${errorPrefix}: ${error.message}`);
    }
}

/**
 * 下载分子图像
 */
function downloadMoleculeImage() {
    try {
        const imgEl = document.getElementById('molecule-img');
        if (!imgEl.src) {
            const noImageMsg = window.i18n ? window.i18n.t('messages.noImage') : '没有可下载的图片';
            alert(noImageMsg);
            return;
        }

        const link = document.createElement('a');
        link.href = imgEl.src;
        link.download = `molecule_prediction_${Date.now()}.png`;
        link.click();

    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.downloadFailed') : '下载失败';
        console.error(errorPrefix, error);
        alert(`${errorPrefix}: ${error.message}`);
    }
}

/**
 * 导出为 CSV 格式（客户端生成）
 *
 * @param {Object} data - 预测数据
 * @param {string} filename - 文件名
 */
function downloadAsCSV(data, filename = 'prediction.csv') {
    try {
        let csv = 'Atom Index,Element,Chemical Shift (ppm)\n';

        data.predictions.forEach(pred => {
            csv += `${pred.atom_index},${pred.element},${pred.ppm.toFixed(2)}\n`;
        });

        csv += '\n';
        csv += `Canonical SMILES,${data.canonical_smiles}\n`;
        csv += `Timestamp,${data.timestamp}\n`;

        downloadFile(csv, filename, 'text/csv');

    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.csvExportFailed') : 'CSV 导出失败';
        console.error(errorPrefix, error);
        alert(`${errorPrefix}: ${error.message}`);
    }
}

/**
 * 导出为 JSON 格式（客户端生成）
 *
 * @param {Object} data - 预测数据
 * @param {string} filename - 文件名
 */
function downloadAsJSON(data, filename = 'prediction.json') {
    try {
        const jsonStr = JSON.stringify(data, null, 2);
        downloadFile(jsonStr, filename, 'application/json');

    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.jsonExportFailed') : 'JSON 导出失败';
        console.error(errorPrefix, error);
        alert(`${errorPrefix}: ${error.message}`);
    }
}

/**
 * 通用文件下载函数
 *
 * @param {string} content - 文件内容
 * @param {string} filename - 文件名
 * @param {string} mimeType - MIME 类型
 */
function downloadFile(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = url;
    link.download = filename;
    link.click();
    URL.revokeObjectURL(url);
}
