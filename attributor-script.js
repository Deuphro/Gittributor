// Attributor Guide - Translation and Interactivity System
// Version 2.0 - With Language Toggle

// Global variables
let currentLanguage = 'fr';
let translations = {};

// Load translations on page load
async function loadTranslations() {
    try {
        const response = await fetch('translations.json');
        if (response.ok) {
            translations = await response.json();
            applyTranslations();
        }
    } catch (error) {
        console.error('Error loading translations:', error);
    }
}

// Apply translations to the page
function applyTranslations() {
    if (!translations || !translations[currentLanguage]) return;
    
    const langData = translations[currentLanguage];
    
    // Translate navigation
    if (langData.nav && translations.meta && translations.meta.languages) {
        translations.meta.languages.forEach(lang => {
            if (langData.nav[lang]) {
                Object.entries(langData.nav[lang]).forEach(([id, text]) => {
                    const navItem = document.querySelector(`a[href="#${id}"]`);
                    if (navItem && navItem.dataset.original) {
                        // Store original text if not already stored
                        if (!navItem.dataset.original) {
                            navItem.dataset.original = navItem.textContent;
                        }
                    }
                });
            }
        });
    }
    
    // Update all translatable elements
    updateTranslatableElements();
    
    // Update language button
    updateLanguageButton();
}

// Update all elements with translation attributes
function updateTranslatableElements() {
    const langData = translations[currentLanguage];
    
    // Process each section
    if (langData.hero) {
        updateHero(langData.hero);
    }
    if (langData.introduction) {
        updateSection('introduction', langData.introduction);
    }
    if (langData.prerequisites) {
        updateSection('prerequisites', langData.prerequisites);
    }
    if (langData.installation) {
        updateSection('installation', langData.installation);
    }
    if (langData.interface) {
        updateSection('interface', langData.interface);
    }
    if (langData.workflow) {
        updateSection('workflow', langData.workflow);
    }
    if (langData.advanced) {
        updateSection('advanced', langData.advanced);
    }
    if (langData.algorithms) {
        updateSection('algorithms', langData.algorithms);
    }
    if (langData.tutorial) {
        updateSection('tutorial', langData.tutorial);
    }
    if (langData.troubleshooting) {
        updateSection('troubleshooting', langData.troubleshooting);
    }
    if (langData.license) {
        updateSection('license', langData.license);
    }
    if (langData.footer) {
        updateFooter(langData.footer);
    }
    if (langData.buttons) {
        updateLanguageButtonText(langData.buttons);
    }
}

// Update hero section
function updateHero(data) {
    const heroSection = document.querySelector('.hero');
    if (heroSection) {
        const h2 = heroSection.querySelector('h2');
        const p = heroSection.querySelector('p');
        
        if (h2 && data.title) {
            h2.textContent = data.title;
        }
        if (p && data.subtitle) {
            p.textContent = data.subtitle;
        }
    }
}

// Update section content
function updateSection(sectionId, data) {
    const section = document.getElementById(sectionId);
    if (!section) return;
    
    // Update section title
    const h2 = section.querySelector('h2');
    if (h2 && data.title) {
        h2.textContent = data.title;
    }
    
    // Section-specific updates
    switch(sectionId) {
        case 'introduction':
            updateIntroductionSection(section, data);
            break;
        case 'prerequisites':
            updatePrerequisitesSection(section, data);
            break;
        case 'installation':
            updateInstallationSection(section, data);
            break;
        case 'interface':
            updateInterfaceSection(section, data);
            break;
        case 'workflow':
            // Workflow section - basic translation support
            break;
        case 'advanced':
            // Advanced section - basic translation support
            break;
        case 'algorithms':
            // Algorithms section - basic translation support
            break;
        case 'tutorial':
            // Tutorial section - basic translation support
            break;
        case 'troubleshooting':
            // Troubleshooting section - basic translation support
            break;
        case 'license':
            // License section - basic translation support
            break;
    }
}

// Update introduction section
function updateIntroductionSection(section, data) {
    const p1 = section.querySelector('.section > p');
    if (p1 && data.p1) {
        p1.textContent = data.p1;
    }
    
    // Update list items
    const list = section.querySelector('ul');
    if (list && data.list) {
        const items = list.querySelectorAll('li');
        data.list.forEach((text, index) => {
            if (items[index]) {
                items[index].textContent = text;
            }
        });
    }
    
    // Update info box
    const infoBox = section.querySelector('.info-box');
    if (infoBox && data.domain_title) {
        const h4 = infoBox.querySelector('h4');
        const p = infoBox.querySelector('p');
        if (h4) h4.textContent = data.domain_title;
        if (p) p.textContent = data.domain_text;
    }
}

// Update prerequisites section
function updatePrerequisitesSection(section, data) {
    // Update steps
    const steps = section.querySelectorAll('.step');
    if (steps && data.steps) {
        Object.entries(data.steps).forEach(([key, stepData]) => {
            const step = steps[parseInt(key) - 1];
            if (step) {
                const h4 = step.querySelector('h4');
                const p = step.querySelector('p');
                if (h4) h4.textContent = stepData.title;
                if (p) p.textContent = stepData.desc;
            }
        });
    }
    
    // Update info box
    const infoBox = section.querySelector('.info-box');
    if (infoBox && data.info_title) {
        const h4 = infoBox.querySelector('h4');
        const p = infoBox.querySelector('p');
        if (h4) h4.textContent = data.info_title;
        if (p) p.textContent = data.info_text;
    }
}

// Update installation section
function updateInstallationSection(section, data) {
    // Update files table
    const table = section.querySelector('table');
    if (table && data.files_table) {
        // This would require more complex DOM manipulation
        // For now, we'll skip table translation for simplicity
    }
    
    // Update steps
    const steps = section.querySelectorAll('.step');
    if (steps && data.steps) {
        Object.entries(data.steps).forEach(([key, stepData]) => {
            const step = steps[parseInt(key) - 1];
            if (step) {
                const h4 = step.querySelector('h4');
                const p = step.querySelector('p');
                if (h4) h4.textContent = stepData.title;
                if (p) p.textContent = stepData.desc;
            }
        });
    }
    
    // Update text sections
    const firstLaunch = section.querySelector('.section > h3');
    if (firstLaunch && data.first_launch) {
        firstLaunch.textContent = data.first_launch;
    }
    
    const windowsText = section.querySelector('.section > p');
    if (windowsText && data.windows_text) {
        windowsText.textContent = data.windows_text;
    }
    
    // Update warning box
    const warningBox = section.querySelector('.warning-box');
    if (warningBox && data.warning_title) {
        const h4 = warningBox.querySelector('h4');
        const p = warningBox.querySelector('p');
        if (h4) h4.textContent = data.warning_title;
        if (p) p.textContent = data.warning_text;
    }
}

// Update interface section
function updateInterfaceSection(section, data) {
    // Update screenshot description
    const screenshotTitle = section.querySelector('.screenshot-section > h3');
    if (screenshotTitle && data.screenshot_title) {
        screenshotTitle.textContent = data.screenshot_title;
    }
    
    const screenshotDesc = section.querySelector('.screenshot-section > p');
    if (screenshotDesc && data.screenshot_desc) {
        screenshotDesc.textContent = data.screenshot_desc;
    }
    
    // Update main windows title
    const mainWindowsTitle = section.querySelector('.section > h3:nth-of-type(2)');
    if (mainWindowsTitle && data.main_windows) {
        mainWindowsTitle.textContent = data.main_windows;
    }
    
    // Update interface zones
    const interfaceZones = section.querySelectorAll('.interface-zone');
    if (interfaceZones && data.panel_desc) {
        const descriptions = [data.panel_desc, data.elaborator_desc, data.agregator_desc, data.dmvm_desc];
        interfaceZones.forEach((zone, index) => {
            const p = zone.querySelector('p');
            if (p && descriptions[index]) {
                p.textContent = descriptions[index];
            }
        });
    }
}

// Update footer section
function updateFooter(data) {
    const footer = document.querySelector('footer');
    if (!footer) return;
    
    // Update copyright
    const copyright = footer.querySelector('.copyright');
    if (copyright && data.copyright) {
        const p = copyright.querySelector('p');
        if (p) p.textContent = data.copyright;
    }
    
    // Update footer left
    const footerLeft = footer.querySelector('.footer-left');
    if (footerLeft && data.left_title) {
        const h3 = footerLeft.querySelector('h3');
        if (h3) h3.textContent = data.left_title;
    }
    
    // Update footer links
    const footerLinks = footer.querySelector('.footer-links');
    if (footerLinks && data.links) {
        const links = footerLinks.querySelectorAll('a');
        data.links.forEach((text, index) => {
            if (links[index]) {
                links[index].textContent = text;
            }
        });
    }
}

// Update language button text
function updateLanguageButtonText(data) {
    const langButton = document.getElementById('language-toggle');
    if (langButton && data[currentLanguage]) {
        langButton.textContent = data[currentLanguage].toggle;
    }
}

// Update language button appearance
function updateLanguageButton() {
    const langButton = document.getElementById('language-toggle');
    if (langButton && translations.buttons && translations.buttons[currentLanguage]) {
        langButton.textContent = translations.buttons[currentLanguage].toggle;
    }
}

// Toggle language
function toggleLanguage() {
    currentLanguage = currentLanguage === 'fr' ? 'en' : 'fr';
    localStorage.setItem('attributor-lang', currentLanguage);
    applyTranslations();
    
    // Update document language attribute
    document.documentElement.lang = currentLanguage;
}

// Smooth scrolling for navigation
document.addEventListener('DOMContentLoaded', () => {
    // Load saved language preference
    const savedLang = localStorage.getItem('attributor-lang');
    if (savedLang && ['fr', 'en'].includes(savedLang)) {
        currentLanguage = savedLang;
    }
    
    // Load translations
    loadTranslations();
    
    // Set up navigation
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                target.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
                
                // Update active nav item
                document.querySelectorAll('.nav-item').forEach(item => {
                    item.classList.remove('active');
                });
                this.classList.add('active');
            }
        });
    });
    
    // Set up language toggle button
    const langButton = document.getElementById('language-toggle');
    if (langButton) {
        langButton.addEventListener('click', toggleLanguage);
    }
    
    // Initialize active nav item on page load
    const firstNavItem = document.querySelector('.nav-item');
    if (firstNavItem) {
        firstNavItem.classList.add('active');
    }
    
    // Add 'no-print' class to elements that shouldn't print
    document.querySelectorAll('nav, .no-print').forEach(el => {
        el.classList.add('no-print');
    });
});

// Highlight nav item on scroll
window.addEventListener('scroll', () => {
    const sections = document.querySelectorAll('.section');
    const navItems = document.querySelectorAll('.nav-item');
    let current = '';

    sections.forEach(section => {
        const sectionTop = section.offsetTop;
        const sectionHeight = section.clientHeight;
        if (window.scrollY >= (sectionTop - 100)) {
            current = section.getAttribute('id');
        }
    });

    navItems.forEach(item => {
        item.classList.remove('active');
        if (item.getAttribute('href') === `#${current}`) {
            item.classList.add('active');
        }
    });
});

// Scroll to top button
function createScrollToTopButton() {
    const button = document.createElement('button');
    button.innerHTML = '↑ Top';
    button.style.cssText = `
        position: fixed;
        bottom: 20px;
        right: 20px;
        padding: 10px 15px;
        background: var(--ipag-light-blue, #3b82f6);
        color: white;
        border: none;
        border-radius: 50px;
        cursor: pointer;
        box-shadow: 0 2px 10px rgba(0,0,0,0.2);
        display: none;
        z-index: 1000;
        font-weight: bold;
    `;
    button.addEventListener('click', () => {
        window.scrollTo({ top: 0, behavior: 'smooth' });
    });
    document.body.appendChild(button);
    
    window.addEventListener('scroll', () => {
        if (window.scrollY > 300) {
            button.style.display = 'block';
        } else {
            button.style.display = 'none';
        }
    });
}

// Call the function to create scroll to top button
document.addEventListener('DOMContentLoaded', createScrollToTopButton);

// Keyboard support
document.addEventListener('keydown', (e) => {
    if (e.key === 'Escape') {
        window.scrollTo({ top: 0, behavior: 'smooth' });
    }
});
