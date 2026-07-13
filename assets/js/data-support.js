// Data support page: the 30-cytokine effect-vs-selectivity showcase.

// one small panel: all genes at (selectivity, logFC on this cytokine); significant hits highlighted
function cytokinePanel(p, ext) {
  const W = 208, H = 176, padL = 30, padR = 8, padT = 20, padB = 26;
  const M = ext.lfc_cap, sMax = ext.sel_max * 1.03;
  const sx = x => padL + x / sMax * (W - padL - padR);
  const sy = y => H - padB - (Math.max(-M, Math.min(M, y)) + M) / (2 * M) * (H - padT - padB);
  const y0 = sy(0);
  const bg = p.points.filter(d => !d.sig).map(d =>
    `<circle cx="${sx(d.sel).toFixed(1)}" cy="${sy(d.lfc).toFixed(1)}" r="1.8" style="fill:var(--muted)" opacity="0.28"/>`).join('');
  const hi = p.points.filter(d => d.sig).map(d => {
    const up = d.lfc >= 0;
    return `<circle cx="${sx(d.sel).toFixed(1)}" cy="${sy(d.lfc).toFixed(1)}" r="3.4"
      style="fill:${up ? 'var(--clay)' : 'var(--eff-down)'};stroke:var(--panel)" stroke-width="0.7">
      <title>${d.gene} — selectivity ${d.sel.toFixed(2)}, log₂FC ${d.lfc >= 0 ? '+' : ''}${d.lfc.toFixed(2)}</title></circle>`;
  }).join('');
  return `<div class="cyt-panel">
    <svg viewBox="0 0 ${W} ${H}" width="100%">
      <text x="${padL}" y="13" class="cyt-title">${p.cytokine}</text>
      <text x="${W - padR}" y="13" text-anchor="end" class="cyt-n">${p.n_sig} hit${p.n_sig === 1 ? '' : 's'}</text>
      <line x1="${padL}" y1="${(H - padB).toFixed(1)}" x2="${(W - padR).toFixed(1)}" y2="${(H - padB).toFixed(1)}" style="stroke:var(--line)"/>
      <line x1="${padL}" y1="${padT}" x2="${padL}" y2="${(H - padB).toFixed(1)}" style="stroke:var(--line)"/>
      <line x1="${padL}" y1="${y0.toFixed(1)}" x2="${(W - padR).toFixed(1)}" y2="${y0.toFixed(1)}" style="stroke:var(--line)" stroke-dasharray="3 3"/>
      <text x="${(padL - 4).toFixed(1)}" y="${(padT + 4).toFixed(1)}" text-anchor="end" class="cyt-tick">+${M}</text>
      <text x="${(padL - 4).toFixed(1)}" y="${(H - padB).toFixed(1)}" text-anchor="end" class="cyt-tick">−${M}</text>
      ${bg}${hi}
    </svg></div>`;
}

function renderShowcase(box, data) {
  const ext = { lfc_cap: data.lfc_cap, sel_max: data.sel_max };
  const panels = data.panels.map(p => cytokinePanel(p, ext)).join('');
  box.innerHTML = `
    <h1 class="showcase-title">Data support — context-specific regulation, cytokine by cytokine</h1>
    <p class="showcase-sub">Each panel is one of the 30 cytokines in the module. Every point is one of the ${data.n_genes} regulators, placed by its <b>selectivity index</b> (x, →) and its <b>log₂ fold-change on that cytokine</b> in stimulated cells (y, fixed at ±${data.lfc_cap}). Coloured points are the gene×cytokine pairs that pass the robust-hit gate (FDR&lt;0.05); grey points are non-significant. Panels are ordered by number of robust hits.</p>
    <div class="showcase-legend">
      <span class="sc-item"><span class="sc-dot" style="background:var(--clay)"></span>Significant · induced by knockdown (log₂FC &gt; 0)</span>
      <span class="sc-item"><span class="sc-dot" style="background:var(--eff-down)"></span>Significant · suppressed by knockdown (log₂FC &lt; 0)</span>
      <span class="sc-item"><span class="sc-dot" style="background:var(--muted);opacity:.4"></span>Not significant</span>
    </div>
    <div class="showcase-grid">${panels}</div>`;
}

Promise.all([
  fetch('data/cytokine_panels.json').then(r => r.ok ? r.json() : null),
  fetch('data/index.json').then(r => r.json()).catch(() => ({})),
]).then(([panels, idx]) => {
  document.title = 'Data support — Target Browser';
  if (panels) renderShowcase(document.getElementById('showcase'), panels);
  else document.getElementById('showcase').innerHTML =
    '<div class="loading">Failed to load data — serve over HTTP (see README).</div>';
  if (window.renderFooter && idx.site) window.renderFooter(idx.site);
}).catch(() => {
  document.getElementById('showcase').innerHTML =
    '<div class="loading">Failed to load data — serve over HTTP (see README).</div>';
});
