/**
 * 主 JavaScript 文件 - 应用初始化和全局功能
 */

console.log('应用启动: 硼核磁预测系统');

/**
 * 应用配置
 */
const APP_CONFIG = {
    apiUrl: '/api',
    defaultSolvent: 'CDCl3',
    autoSaveSolvent: true,
    retryAttempts: 3,
    retryDelay: 1000 // ms
};

/**
 * 应用状态
 */
const APP_STATE = {
    isLoading: false,
    currentPrediction: null
};

/**
 * 初始化应用
 */
function initializeApp() {
    console.log('初始化应用...');

    // 恢复保存的溶剂选择
    restoreSolventSelection();

    // 绑定全局事件
    bindGlobalEvents();

    // 初始化 localStorage
    initializeLocalStorage();

    console.log('✓ 应用初始化完成');
}

/**
 * 恢复保存的溶剂选择
 */
function restoreSolventSelection() {
    if (!APP_CONFIG.autoSaveSolvent) return;

    const solventSelect = document.getElementById('solvent-select');
    if (!solventSelect) return;

    const savedSolvent = localStorage.getItem('selectedSolvent');
    if (savedSolvent && solventSelect.options.namedItem(savedSolvent)) {
        solventSelect.value = savedSolvent;
    }

    // 保存溶剂选择
    solventSelect.addEventListener('change', () => {
        localStorage.setItem('selectedSolvent', solventSelect.value);
    });
}

/**
 * 初始化 localStorage
 */
function initializeLocalStorage() {
    // 检查浏览器是否支持 localStorage
    if (typeof localStorage === 'undefined') {
        console.warn('浏览器不支持 localStorage');
        return;
    }

    console.log('✓ localStorage 已初始化');
}

/**
 * 绑定全局事件
 */
function bindGlobalEvents() {
    // 处理网络连接状态
    window.addEventListener('online', () => {
        console.log('✓ 网络已连接');
        const message = window.i18n ? window.i18n.t('messages.networkOnline') : '网络已连接';
        showNotification(message, 'success', 2000);
    });

    window.addEventListener('offline', () => {
        console.log('✗ 网络已断开');
        const message = window.i18n ? window.i18n.t('messages.networkOffline') : '网络已断开';
        showNotification(message, 'error', 3000);
    });

    // 处理页面可见性
    document.addEventListener('visibilitychange', () => {
        if (document.visibilityState === 'visible') {
            console.log('页面重新可见');
        } else {
            console.log('页面隐藏');
        }
    });
}

/**
 * 显示通知消息
 *
 * @param {string} message - 消息内容
 * @param {string} type - 消息类型: 'success', 'error', 'warning', 'info'
 * @param {number} duration - 显示时长 (ms), 0 = 不自动关闭
 */
function showNotification(message, type = 'info', duration = 3000) {
    // 创建通知容器（如果不存在）
    let notificationContainer = document.getElementById('notification-container');
    if (!notificationContainer) {
        notificationContainer = document.createElement('div');
        notificationContainer.id = 'notification-container';
        notificationContainer.style.cssText = `
            position: fixed;
            top: 20px;
            right: 20px;
            z-index: 9999;
            max-width: 400px;
        `;
        document.body.appendChild(notificationContainer);
    }

    // 创建通知元素
    const notification = document.createElement('div');
    notification.className = `notification notification-${type}`;
    notification.textContent = message;
    notification.style.cssText = `
        padding: 12px 16px;
        margin-bottom: 10px;
        border-radius: 4px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.15);
        animation: slideIn 0.3s ease;
        font-weight: 500;
    `;

    // 根据类型设置样式
    const colors = {
        success: { bg: '#48bb78', color: '#fff' },
        error: { bg: '#f56565', color: '#fff' },
        warning: { bg: '#ed8936', color: '#fff' },
        info: { bg: '#4299e1', color: '#fff' }
    };

    const style = colors[type] || colors.info;
    notification.style.backgroundColor = style.bg;
    notification.style.color = style.color;

    // 添加到容器
    notificationContainer.appendChild(notification);

    // 自动关闭
    if (duration > 0) {
        setTimeout(() => {
            notification.style.animation = 'slideOut 0.3s ease';
            setTimeout(() => notification.remove(), 300);
        }, duration);
    }

    return notification;
}

/**
 * 添加动画样式
 */
function addAnimationStyles() {
    const styleSheet = document.createElement('style');
    styleSheet.textContent = `
        @keyframes slideIn {
            from {
                transform: translateX(400px);
                opacity: 0;
            }
            to {
                transform: translateX(0);
                opacity: 1;
            }
        }

        @keyframes slideOut {
            from {
                transform: translateX(0);
                opacity: 1;
            }
            to {
                transform: translateX(400px);
                opacity: 0;
            }
        }
    `;
    document.head.appendChild(styleSheet);
}

/**
 * 检查网络连接
 *
 * @returns {boolean} 是否有网络连接
 */
function isOnline() {
    return navigator.onLine;
}

/**
 * 获取 API URL
 *
 * @param {string} endpoint - API 端点
 * @returns {string} 完整 API URL
 */
function getApiUrl(endpoint) {
    return `${APP_CONFIG.apiUrl}${endpoint}`;
}

/**
 * 带重试的 Fetch 请求
 *
 * @param {string} url - URL
 * @param {Object} options - 选项
 * @param {number} attempt - 当前尝试次数
 * @returns {Promise<Response>}
 */
async function fetchWithRetry(url, options = {}, attempt = 1) {
    try {
        const response = await fetch(url, options);
        return response;
    } catch (error) {
        if (attempt < APP_CONFIG.retryAttempts) {
            console.warn(`请求失败，${APP_CONFIG.retryDelay}ms 后重试 (${attempt}/${APP_CONFIG.retryAttempts})...`);
            await new Promise(resolve => setTimeout(resolve, APP_CONFIG.retryDelay));
            return fetchWithRetry(url, options, attempt + 1);
        }
        throw error;
    }
}

/**
 * 格式化日期
 *
 * @param {string|Date} date - 日期
 * @param {string} format - 格式字符串
 * @returns {string} 格式化的日期
 */
function formatDate(date, format = 'YYYY-MM-DD HH:mm:ss') {
    const d = new Date(date);
    const year = d.getFullYear();
    const month = String(d.getMonth() + 1).padStart(2, '0');
    const day = String(d.getDate()).padStart(2, '0');
    const hours = String(d.getHours()).padStart(2, '0');
    const minutes = String(d.getMinutes()).padStart(2, '0');
    const seconds = String(d.getSeconds()).padStart(2, '0');

    return format
        .replace('YYYY', year)
        .replace('MM', month)
        .replace('DD', day)
        .replace('HH', hours)
        .replace('mm', minutes)
        .replace('ss', seconds);
}

/**
 * 防抖函数
 *
 * @param {Function} func - 要防抖的函数
 * @param {number} wait - 等待时间 (ms)
 * @returns {Function} 防抖后的函数
 */
function debounce(func, wait = 300) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
}

/**
 * 节流函数
 *
 * @param {Function} func - 要节流的函数
 * @param {number} limit - 限制时间 (ms)
 * @returns {Function} 节流后的函数
 */
function throttle(func, limit = 300) {
    let inThrottle;
    return function(...args) {
        if (!inThrottle) {
            func.apply(this, args);
            inThrottle = true;
            setTimeout(() => inThrottle = false, limit);
        }
    };
}

/**
 * 页面加载完成后初始化
 */
document.addEventListener('DOMContentLoaded', () => {
    addAnimationStyles();
    initializeApp();

    // 开发环境日志
    if (window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1') {
        console.log('开发模式: 已启用详细日志');
    }
});

/**
 * 页面卸载前清理
 */
window.addEventListener('beforeunload', () => {
    console.log('页面即将卸载');
});

// 导出全局函数供其他脚本使用
window.showNotification = showNotification;
window.fetchWithRetry = fetchWithRetry;
window.formatDate = formatDate;
window.debounce = debounce;
window.throttle = throttle;
window.isOnline = isOnline;
window.getApiUrl = getApiUrl;
