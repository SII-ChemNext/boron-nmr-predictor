/**
 * Lightweight i18n Translation Engine
 * Handles loading, caching, and applying translations
 */

class I18n {
    constructor() {
        this.currentLang = this.detectLanguage();
        this.translations = {};
        this.fallbackLang = 'zh-CN';
    }

    /**
     * Detect user's preferred language
     */
    detectLanguage() {
        // 1. Check localStorage
        const saved = localStorage.getItem('preferredLanguage');
        if (saved) {
            console.log(`使用保存的语言偏好: ${saved}`);
            return saved;
        }

        // 2. Check browser language
        const browserLang = navigator.language || navigator.userLanguage;
        console.log(`浏览器语言: ${browserLang}`);

        if (browserLang.startsWith('zh')) {
            return 'zh-CN';
        }
        if (browserLang.startsWith('en')) {
            return 'en-US';
        }

        // 3. Default to Chinese
        return 'zh-CN';
    }

    /**
     * Load translation file
     */
    async loadTranslations(lang) {
        if (this.translations[lang]) {
            console.log(`翻译已缓存: ${lang}`);
            return this.translations[lang];
        }

        try {
            console.log(`加载翻译文件: ${lang}.json`);
            const response = await fetch((window.BASE_PATH || '/') + `static/i18n/${lang}.json`);

            if (!response.ok) {
                throw new Error(`HTTP ${response.status}`);
            }

            const data = await response.json();
            this.translations[lang] = data;
            console.log(`✓ 翻译文件加载成功: ${lang}`);
            return data;

        } catch (error) {
            console.error(`加载翻译文件失败 (${lang}):`, error);

            // Try fallback language
            if (lang !== this.fallbackLang && !this.translations[this.fallbackLang]) {
                console.log(`尝试加载备用语言: ${this.fallbackLang}`);
                return this.loadTranslations(this.fallbackLang);
            }

            return {};
        }
    }

    /**
     * Get translation by key path
     * Example: i18n.t('app.title')
     */
    t(keyPath, defaultValue = '') {
        const translation = this.translations[this.currentLang];
        if (!translation) {
            console.warn(`翻译未加载: ${this.currentLang}`);
            return defaultValue || keyPath;
        }

        const keys = keyPath.split('.');
        let value = translation;

        for (const key of keys) {
            if (value && typeof value === 'object' && key in value) {
                value = value[key];
            } else {
                console.warn(`翻译键不存在: ${keyPath}`);
                return defaultValue || keyPath;
            }
        }

        return value;
    }

    /**
     * Change current language
     */
    async changeLanguage(lang) {
        if (lang === this.currentLang) {
            console.log(`语言已经是: ${lang}`);
            return;
        }

        console.log(`切换语言: ${this.currentLang} -> ${lang}`);

        await this.loadTranslations(lang);
        this.currentLang = lang;
        localStorage.setItem('preferredLanguage', lang);

        // Update all translatable elements
        this.updateDOM();

        // Dispatch event for custom handlers
        window.dispatchEvent(new CustomEvent('languageChanged', {
            detail: { language: lang }
        }));

        console.log(`✓ 语言切换完成: ${lang}`);
    }

    /**
     * Update DOM with current language
     */
    updateDOM() {
        console.log('更新DOM元素...');
        let updatedCount = 0;

        // Update elements with data-i18n attribute
        document.querySelectorAll('[data-i18n]').forEach(el => {
            const key = el.getAttribute('data-i18n');
            const translation = this.t(key);

            if (translation && translation !== key) {
                // Check if we should update text or placeholder
                if (el.tagName === 'INPUT' || el.tagName === 'TEXTAREA') {
                    if (el.hasAttribute('placeholder')) {
                        el.placeholder = translation;
                    } else {
                        el.value = translation;
                    }
                } else if (el.tagName === 'IMG') {
                    el.alt = translation;
                } else {
                    el.textContent = translation;
                }
                updatedCount++;
            }
        });

        // Update elements with data-i18n-html (for HTML content)
        document.querySelectorAll('[data-i18n-html]').forEach(el => {
            const key = el.getAttribute('data-i18n-html');
            const translation = this.t(key);
            if (translation && translation !== key) {
                el.innerHTML = translation;
                updatedCount++;
            }
        });

        // Update elements with data-i18n-title attribute (for title tooltips)
        document.querySelectorAll('[data-i18n-title]').forEach(el => {
            const key = el.getAttribute('data-i18n-title');
            const translation = this.t(key);
            if (translation && translation !== key) {
                el.title = translation;
                updatedCount++;
            }
        });

        // Update document title
        const titleKey = document.documentElement.getAttribute('data-i18n-title');
        if (titleKey) {
            const titleTranslation = this.t(titleKey);
            if (titleTranslation && titleTranslation !== titleKey) {
                document.title = titleTranslation;
                updatedCount++;
            }
        }

        // Update language attribute
        document.documentElement.lang = this.currentLang === 'zh-CN' ? 'zh-CN' : 'en-US';

        console.log(`✓ 更新了 ${updatedCount} 个DOM元素`);
    }

    /**
     * Initialize i18n system
     */
    async init() {
        console.log('初始化 i18n 系统...');
        console.log(`当前语言: ${this.currentLang}`);

        await this.loadTranslations(this.currentLang);
        this.updateDOM();

        console.log('✓ i18n 初始化完成');
    }

    /**
     * Get current language code
     */
    getCurrentLanguage() {
        return this.currentLang;
    }

    /**
     * Get display name for language
     */
    getLanguageName(lang = this.currentLang) {
        const names = {
            'zh-CN': '中文',
            'en-US': 'English'
        };
        return names[lang] || lang;
    }
}

// Create global instance
window.i18n = new I18n();

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    window.i18n.init();
});
