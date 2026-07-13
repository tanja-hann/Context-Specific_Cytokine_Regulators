// Theme: light default, persisted, applied before paint (set inline in <head>).
(function () {
  const KEY = 'tb-theme';
  function apply(t) {
    document.documentElement.setAttribute('data-theme', t);
    const btn = document.getElementById('theme-toggle');
    if (btn) btn.innerHTML = (t === 'dark')
      ? '<span>☀︎</span> Light'
      : '<span>☾</span> Dark';
    // let charts re-read CSS variables
    window.dispatchEvent(new CustomEvent('themechange', { detail: t }));
  }
  window.__getTheme = () => document.documentElement.getAttribute('data-theme') || 'light';
  window.__initTheme = function () {
    const saved = localStorage.getItem(KEY) || 'light';   // default light
    apply(saved);
    const btn = document.getElementById('theme-toggle');
    if (btn) btn.addEventListener('click', () => {
      const next = window.__getTheme() === 'dark' ? 'light' : 'dark';
      localStorage.setItem(KEY, next);
      apply(next);
    });
  };
  // read a CSS custom property from :root (theme-aware)
  window.cssVar = name => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  // shared site footer (attribution + citation), rendered from index.json's site block
  window.renderFooter = function (site) {
    const el = document.getElementById('site-footer');
    if (!el || !site) return;
    const esc = s => (s || '').replace(/[&<>"]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
    let html = '';
    if (site.attribution) html += `<div class="foot-attr">${esc(site.attribution)}</div>`;
    if (site.footer) html += `<div class="foot-note">${esc(site.footer)}</div>`;
    if (site.repo_url) html += `<div class="foot-repo"><a href="${esc(site.repo_url)}" target="_blank" rel="noopener">⟩ Source code on GitHub</a></div>`;
    if (site.citation && site.citation.text) {
      const t = esc(site.citation.text);
      const cited = site.citation.url
        ? `${t} <a href="${esc(site.citation.url)}" target="_blank" rel="noopener">${esc(site.citation.url)}</a>`
        : t;
      html += `<div class="foot-cite"><span class="foot-cite-label">Cite</span>${cited}</div>`;
    }
    el.innerHTML = html;
  };
})();
