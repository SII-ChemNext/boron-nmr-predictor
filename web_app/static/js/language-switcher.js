/**
 * Language Switcher Component
 * Handles language toggle button and switching logic
 */

class LanguageSwitcher {
    constructor() {
        this.button = null;
        this.isToggling = false;
    }

    /**
     * Initialize the language switcher
     */
    init() {
        // Wait for i18n to be ready
        if (!window.i18n) {
            console.warn('i18n not available, retrying in 100ms...');
            setTimeout(() => this.init(), 100);
            return;
        }

        this.createButton();
        this.bindEvents();
        this.updateButtonText();

        console.log('✓ Language switcher initialized');
    }

    /**
     * Create language toggle button
     */
    createButton() {
        // Find the nav element
        const nav = document.querySelector('.nav');
        if (!nav) {
            console.warn('Nav element not found, cannot create language switcher');
            return;
        }

        // Create button
        this.button = document.createElement('button');
        this.button.id = 'language-toggle';
        this.button.className = 'btn-language';
        this.button.innerHTML = `
            <span class="lang-icon">🌐</span>
            <span class="lang-text">中/EN</span>
        `;
        this.button.setAttribute('title', 'Switch Language / 切换语言');

        // Append to nav
        nav.appendChild(this.button);
    }

    /**
     * Bind event listeners
     */
    bindEvents() {
        if (!this.button) return;

        this.button.addEventListener('click', async () => {
            if (this.isToggling) return;

            this.isToggling = true;
            this.button.disabled = true;

            try {
                const currentLang = window.i18n.getCurrentLanguage();
                const newLang = currentLang === 'zh-CN' ? 'en-US' : 'zh-CN';

                console.log(`Switching language: ${currentLang} -> ${newLang}`);

                await window.i18n.changeLanguage(newLang);
                this.updateButtonText();

                // Show notification
                const message = newLang === 'zh-CN'
                    ? window.i18n.t('messages.languageSwitched')
                    : window.i18n.t('messages.languageSwitched');

                if (window.showNotification) {
                    window.showNotification(message, 'success', 2000);
                }

            } catch (error) {
                console.error('Language switch failed:', error);

                const errorMsg = window.i18n ? window.i18n.t('messages.languageSwitchFailed') : 'Language switch failed / 语言切换失败';

                if (window.showNotification) {
                    window.showNotification(errorMsg, 'error', 3000);
                }
            } finally {
                this.isToggling = false;
                this.button.disabled = false;
            }
        });

        // Listen for language change events
        window.addEventListener('languageChanged', () => {
            this.updateButtonText();
        });
    }

    /**
     * Update button text based on current language
     */
    updateButtonText() {
        if (!this.button) return;

        const currentLang = window.i18n.getCurrentLanguage();
        const langText = this.button.querySelector('.lang-text');

        if (langText) {
            // Show opposite language as the target
            langText.textContent = currentLang === 'zh-CN' ? 'EN' : '中';
        }
    }
}

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', () => {
    const switcher = new LanguageSwitcher();
    switcher.init();
});
