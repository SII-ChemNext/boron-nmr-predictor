/**
 * Ketcher 分子编辑器集成
 * 直接访问 iframe 中的 Ketcher API
 */

let ketcherInstance = null;
let ketcherIframe = null;

/**
 * 初始化 Ketcher 编辑器
 */
async function initializeKetcher() {
    try {
        console.log('初始化 Ketcher 编辑器...');

        ketcherIframe = document.getElementById('ketcher-iframe');
        if (!ketcherIframe) {
            console.error('找不到 Ketcher iframe');
            createFallbackInput();
            return;
        }

        // 等待 iframe 加载完成
        console.log('等待 iframe 加载...');
        await waitForIframeReady(15000);

        // 检查 Ketcher API 是否可用
        if (!ketcherIframe.contentWindow || !ketcherIframe.contentWindow.ketcher) {
            console.warn('Ketcher API 不可用，使用备选输入框');
            createFallbackInput();
            return;
        }

        console.log('✓ Ketcher iframe 已就绪');

        // 设置 ketcherInstance
        ketcherInstance = ketcherIframe.contentWindow.ketcher;
        console.log('✓ Ketcher 编辑器已初始化成功');

        // 检查是否需要加载保存的分子
        const savedMolecule = localStorage.getItem('loadMolecule');
        if (savedMolecule) {
            try {
                await setSmilesInKetcher(savedMolecule);
                localStorage.removeItem('loadMolecule');
                console.log('✓ 已加载保存的分子结构');
            } catch (error) {
                console.warn('加载保存的分子结构失败:', error);
            }
        }

    } catch (error) {
        console.error('Ketcher 初始化失败:', error);
        createFallbackInput();
    }
}

/**
 * 等待 iframe 加载完成
 */
function waitForIframeReady(timeout = 15000) {
    return new Promise((resolve) => {
        let elapsed = 0;

        const checkInterval = setInterval(() => {
            try {
                // 检查 Ketcher 对象是否已加载
                if (ketcherIframe &&
                    ketcherIframe.contentWindow &&
                    ketcherIframe.contentWindow.ketcher) {
                    console.log('✓ Ketcher 库已加载');
                    clearInterval(checkInterval);
                    resolve(true);
                    return;
                }
            } catch (e) {
                // iframe 还未完全加载
            }

            elapsed += 100;
            if (elapsed > timeout) {
                clearInterval(checkInterval);
                console.warn('Ketcher 加载超时');
                resolve(false);
            }
        }, 100);
    });
}

/**
 * 创建备选的 SMILES 输入框
 */
function createFallbackInput() {
    const wrapper = document.getElementById('ketcher-wrapper');
    if (wrapper) {
        const fallbackText = window.i18n ? window.i18n.t('editor.ketcherFallback') : '⚠️ Ketcher 编辑器加载失败，请直接输入 SMILES';
        const placeholderText = window.i18n ? window.i18n.t('editor.smilesPlaceholder') : '输入分子 SMILES，例如：OB(O)c1ccccc1';

        wrapper.innerHTML = `
            <div style="padding: 20px; background: white; height: 100%; display: flex; flex-direction: column; justify-content: center; align-items: center;">
                <p style="color: #666; margin-bottom: 15px; font-size: 14px;">
                    ${fallbackText}
                </p>
                <input type="text" id="smiles-fallback-input"
                       placeholder="${placeholderText}"
                       style="width: 100%; padding: 10px; border: 1px solid #ccc; border-radius: 4px; font-family: monospace; font-size: 14px;">
            </div>
        `;
    }
}

/**
 * 从 Ketcher 获取 SMILES
 */
async function getSmilesFromKetcher() {
    try {
        // 检查是否使用备选输入框
        const fallbackInput = document.getElementById('smiles-fallback-input');
        if (fallbackInput) {
            const smiles = fallbackInput.value.trim();
            if (!smiles) {
                const errorMsg = window.i18n ? window.i18n.t('messages.inputSmiles') : '请输入 SMILES 字符串';
                throw new Error(errorMsg);
            }
            return smiles;
        }

        // 使用 Ketcher API
        if (!ketcherInstance) {
            const errorMsg = window.i18n ? window.i18n.t('messages.ketcherNotInitialized') : 'Ketcher 编辑器未初始化';
            throw new Error(errorMsg);
        }

        console.log('从 Ketcher 获取 SMILES...');

        // 调用 Ketcher 的 getSmiles 方法
        // 这个方法返回一个 Promise
        let smiles = await ketcherInstance.getSmiles();

        if (!smiles || smiles.trim() === '') {
            const errorMsg = window.i18n ? window.i18n.t('messages.drawMolecule') : '请先在编辑器中绘制分子结构';
            throw new Error(errorMsg);
        }

        console.log('✓ 获得 SMILES:', smiles);
        return smiles.trim();

    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.getSmilesError') : '获取 SMILES 失败';
        console.error(errorPrefix, error);
        throw new Error(`${errorPrefix}: ${error.message}`);
    }
}

/**
 * 设置 Ketcher 中的分子
 */
async function setSmilesInKetcher(smiles) {
    try {
        // 检查是否使用备选输入框
        const fallbackInput = document.getElementById('smiles-fallback-input');
        if (fallbackInput) {
            fallbackInput.value = smiles;
            return;
        }

        // 使用 Ketcher API
        if (!ketcherInstance) {
            const errorMsg = window.i18n ? window.i18n.t('messages.ketcherNotInitialized') : 'Ketcher 编辑器未初始化';
            throw new Error(errorMsg);
        }

        console.log('设置分子到 Ketcher:', smiles);

        // 调用 Ketcher 的 setMolecule 方法
        await ketcherInstance.setMolecule(smiles);

        console.log('✓ 分子已设置');

    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.setMoleculeFailed') : '设置分子结构失败';
        console.error(errorPrefix, error);
        throw new Error(`${errorPrefix}: ${error.message}`);
    }
}

/**
 * 清空 Ketcher 画布
 */
function clearKetcher() {
    try {
        // 检查是否使用备选输入框
        const fallbackInput = document.getElementById('smiles-fallback-input');
        if (fallbackInput) {
            fallbackInput.value = '';
            return;
        }

        // 使用 Ketcher API
        if (ketcherInstance) {
            ketcherInstance.setMolecule('');
            console.log('✓ 已清空画布');
        }

    } catch (error) {
        console.error('清空画布失败:', error);
    }
}

/**
 * 加载示例分子
 */
async function loadExampleMolecule() {
    const exampleSmiles = 'OB(O)c1ccccc1'; // 苯硼酸

    try {
        console.log('加载示例分子:', exampleSmiles);
        await setSmilesInKetcher(exampleSmiles);
        console.log('✓ 已加载示例分子');
    } catch (error) {
        const errorPrefix = window.i18n ? window.i18n.t('messages.exampleLoadFailed') : '加载示例分子失败';
        console.error(errorPrefix, error);
        alert(`${errorPrefix}: ${error.message}`);
    }
}

/**
 * 页面加载完成后初始化
 */
document.addEventListener('DOMContentLoaded', async () => {
    console.log('页面已加载，初始化 Ketcher...');

    await initializeKetcher();

    // 绑定按钮事件
    const clearBtn = document.getElementById('clear-btn');
    const exampleBtn = document.getElementById('example-btn');

    if (clearBtn) {
        clearBtn.addEventListener('click', () => {
            clearKetcher();
        });
    }

    if (exampleBtn) {
        exampleBtn.addEventListener('click', loadExampleMolecule);
    }
});
