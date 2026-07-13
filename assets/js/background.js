// Background page: the research write-up (rendered from data/writeup.json).
Promise.all([
  fetch('data/writeup.json').then(r => r.json()),
  fetch('data/index.json').then(r => r.json()).catch(() => ({})),
]).then(([w, idx]) => {
  document.title = (w.title || 'Background') + ' — Target Browser';
  document.getElementById('writeup').innerHTML =
    `<h1 class="writeup-title">${w.title || ''}</h1>` + (w.html || '');
  if (window.renderFooter && idx.site) window.renderFooter(idx.site);
}).catch(() => {
  document.getElementById('writeup').innerHTML =
    '<div class="loading">Failed to load write-up — serve over HTTP (see README).</div>';
});
