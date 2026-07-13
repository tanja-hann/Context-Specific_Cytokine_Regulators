// Landing page: data-driven sortable/filterable table with color-bar scores.
const PALETTE = {
  blue:  ['#5a3a2e', '#d97757'],  /* Claude clay */
  green: ['#4a4230', '#bfa98a'],  /* warm sand */
  amber: ['#5a3f20', '#e0a458'],  /* amber-gold */
  red:   ['#4a2118', '#cc6d52'],  /* deep terracotta */
};
let STATE = { data: null, sortKey: 'symbol', sortDir: 1, q: '', flags: new Set() };

function barHTML(v, cfg) {
  if (v === null || v === undefined)
    return '<div class="barbox"><div class="bar"></div><span class="val na">n/a</span></div>';
  const [c0, c1] = PALETTE[cfg.palette] || PALETTE.blue;
  const frac = Math.max(0, Math.min(1, (v - cfg.min) / (cfg.max - cfg.min)));
  const pct = (frac * 100).toFixed(1);
  return `<div class="barbox" title="${cfg.label}: ${v}">
    <div class="bar"><span style="width:${pct}%;background:linear-gradient(90deg,${c0},${c1})"></span></div>
    <span class="val">${v.toFixed(2)}</span></div>`;
}

// Modality → CSS var (threaded with the rest of the theme). Four distinct modalities.
const MOD_COLOR = {
  'Small molecule':   'var(--blue)',
  'Antibody/biologic':'var(--amber)',
  'ASO':              'var(--green)',
  'siRNA':            'var(--clay)',
  'Undetermined':     'var(--muted)',
};
const MOD_SHAPE = { 'Small molecule':'rect', 'Antibody/biologic':'diamond', 'ASO':'tri', 'siRNA':'circle', 'Undetermined':'circle' };
const QUAD_THRESH = 0.30;  // protein/RNA relative-score cutoff separating the four viability quadrants

function marker(cx, cy, shape, color, r, attrs, child) {
  const a = `style="fill:${color}" ${attrs||''}`;
  const c = child || '';
  if (shape === 'rect')    return `<rect x="${(cx-r).toFixed(1)}" y="${(cy-r).toFixed(1)}" width="${(2*r).toFixed(1)}" height="${(2*r).toFixed(1)}" ${a}>${c}</rect>`;
  if (shape === 'diamond') return `<path d="M ${cx.toFixed(1)} ${(cy-r*1.3).toFixed(1)} L ${(cx+r*1.3).toFixed(1)} ${cy.toFixed(1)} L ${cx.toFixed(1)} ${(cy+r*1.3).toFixed(1)} L ${(cx-r*1.3).toFixed(1)} ${cy.toFixed(1)} Z" ${a}>${c}</path>`;
  if (shape === 'tri')     return `<path d="M ${cx.toFixed(1)} ${(cy-r*1.4).toFixed(1)} L ${(cx+r*1.3).toFixed(1)} ${(cy+r).toFixed(1)} L ${(cx-r*1.3).toFixed(1)} ${(cy+r).toFixed(1)} Z" ${a}>${c}</path>`;
  return `<circle cx="${cx.toFixed(1)}" cy="${cy.toFixed(1)}" r="${r}" ${a}>${c}</circle>`;
}

function renderMatrix(d) {
  const box = document.getElementById('matrix');
  if (!box) return;
  const genes = d.genes.filter(g => g.scores.protein_drug != null && g.scores.rna_drug != null);
  const W = 720, H = 560, padL = 66, padR = 20, padT = 16, padB = 58;
  const x0 = 0, x1 = 1, y0 = 0, y1 = 1;
  const sx = v => padL + (v - x0) / (x1 - x0) * (W - padL - padR);
  const sy = v => H - padB - (v - y0) / (y1 - y0) * (H - padT - padB);
  const ticks = [0, 0.2, 0.4, 0.6, 0.8, 1.0];
  const grid = ticks.map(t => `
    <line x1="${sx(t).toFixed(1)}" y1="${padT}" x2="${sx(t).toFixed(1)}" y2="${(H-padB).toFixed(1)}" style="stroke:var(--line)" stroke-width="0.5" opacity="0.5"/>
    <line x1="${padL}" y1="${sy(t).toFixed(1)}" x2="${(W-padR).toFixed(1)}" y2="${sy(t).toFixed(1)}" style="stroke:var(--line)" stroke-width="0.5" opacity="0.5"/>
    <text x="${sx(t).toFixed(1)}" y="${(H-padB+16).toFixed(1)}" text-anchor="middle" class="mtick">${t.toFixed(1)}</text>
    <text x="${(padL-8).toFixed(1)}" y="${(sy(t)+3.5).toFixed(1)}" text-anchor="end" class="mtick">${t.toFixed(1)}</text>`).join('');
  // four viability quadrants split at QUAD_THRESH on each axis
  const tx = sx(QUAD_THRESH), ty = sy(QUAD_THRESH);
  const qcount = lab => genes.filter(g => {
    const p = g.scores.protein_drug >= QUAD_THRESH, r = g.scores.rna_drug >= QUAD_THRESH;
    return (lab==='Both viable'&&p&&r)||(lab==='Protein-preferred'&&p&&!r)||
           (lab==='RNA-preferred'&&!p&&r)||(lab==='Neither viable'&&!p&&!r);
  }).length;
  const quadShade = `
    <rect x="${sx(0).toFixed(1)}" y="${sy(1).toFixed(1)}" width="${(tx-sx(0)).toFixed(1)}" height="${(ty-sy(1)).toFixed(1)}" fill="var(--green)" opacity="0.05"/>
    <rect x="${tx.toFixed(1)}" y="${sy(1).toFixed(1)}" width="${(sx(1)-tx).toFixed(1)}" height="${(ty-sy(1)).toFixed(1)}" fill="var(--amber)" opacity="0.06"/>
    <rect x="${sx(0).toFixed(1)}" y="${ty.toFixed(1)}" width="${(tx-sx(0)).toFixed(1)}" height="${(sy(0)-ty).toFixed(1)}" fill="var(--muted)" opacity="0.05"/>
    <rect x="${tx.toFixed(1)}" y="${ty.toFixed(1)}" width="${(sx(1)-tx).toFixed(1)}" height="${(sy(0)-ty).toFixed(1)}" fill="var(--blue)" opacity="0.05"/>
    <line x1="${tx.toFixed(1)}" y1="${padT}" x2="${tx.toFixed(1)}" y2="${(H-padB).toFixed(1)}" style="stroke:var(--muted)" stroke-width="1" stroke-dasharray="2 3" opacity="0.55"/>
    <line x1="${padL}" y1="${ty.toFixed(1)}" x2="${(W-padR).toFixed(1)}" y2="${ty.toFixed(1)}" style="stroke:var(--muted)" stroke-width="1" stroke-dasharray="2 3" opacity="0.55"/>`;
  const quadLabels = `
    <text x="${(sx(0)+6).toFixed(1)}" y="${(sy(1)+13).toFixed(1)}" class="mquad">RNA-preferred · ${qcount('RNA-preferred')}</text>
    <text x="${(sx(1)-6).toFixed(1)}" y="${(sy(1)+13).toFixed(1)}" text-anchor="end" class="mquad">Both viable · ${qcount('Both viable')}</text>
    <text x="${(sx(0)+6).toFixed(1)}" y="${(sy(0)-7).toFixed(1)}" class="mquad">Neither viable · ${qcount('Neither viable')}</text>
    <text x="${(sx(1)-6).toFixed(1)}" y="${(sy(0)-7).toFixed(1)}" text-anchor="end" class="mquad">Protein-preferred · ${qcount('Protein-preferred')}</text>`;
  const diag = `<line x1="${sx(0).toFixed(1)}" y1="${sy(0).toFixed(1)}" x2="${sx(1).toFixed(1)}" y2="${sy(1).toFixed(1)}" style="stroke:var(--muted)" stroke-width="1" stroke-dasharray="5 4" opacity="0.5"/>`;
  const pts = genes.map(g => {
    const cx = sx(g.scores.protein_drug), cy = sy(g.scores.rna_drug);
    const color = MOD_COLOR[g.modality] || 'var(--muted)';
    const shape = MOD_SHAPE[g.modality] || 'circle';
    const title = `<title>${g.symbol} — ${g.protein_name||''}\nProtein ${g.scores.protein_drug.toFixed(2)} · RNA ${g.scores.rna_drug.toFixed(2)} · ${g.modality}</title>`;
    return marker(cx, cy, shape, color, 5, `class="mpt" data-sym="${g.symbol}"`, title);
  }).join('');
  const mods = d.modalities.filter(mm => genes.some(g => g.modality === mm));
  // legend rendered as HTML below the plot (keeps it off the point cloud & quadrant labels)
  const legendHTML = '<div class="matrix-legend">' + mods.map(mm => {
    const n = genes.filter(g => g.modality === mm).length;
    const sh = MOD_SHAPE[mm] || 'circle';
    const svgm = `<svg width="14" height="14" viewBox="0 0 14 14">${marker(7,7,sh,MOD_COLOR[mm]||'var(--muted)',5)}</svg>`;
    return `<span class="mleg-item">${svgm}${mm} <b>(${n})</b></span>`;
  }).join('') + '</div>';
  // diagonal label placed ON the diagonal near the middle, rotated to match its slope
  const dmx = sx(0.5), dmy = sy(0.5);
  const dangle = Math.atan2(sy(0.62)-sy(0.5), sx(0.62)-sx(0.5)) * 180 / Math.PI;
  box.innerHTML = `<svg viewBox="0 0 ${W} ${H}" width="100%" class="matrix-svg">
    ${grid}${quadShade}${quadLabels}${diag}
    <text x="${dmx.toFixed(1)}" y="${(dmy-6).toFixed(1)}" text-anchor="middle" class="mdiag" transform="rotate(${dangle.toFixed(1)} ${dmx.toFixed(1)} ${dmy.toFixed(1)})">protein = RNA</text>
    ${pts}
    <text x="${((padL+W-padR)/2).toFixed(1)}" y="${H-16}" text-anchor="middle" class="maxis">Protein-level druggability (relative) →</text>
    <text x="18" y="${((padT+H-padB)/2).toFixed(1)}" text-anchor="middle" transform="rotate(-90 18 ${((padT+H-padB)/2).toFixed(1)})" class="maxis">RNA-level druggability (relative) →</text>
  </svg>${legendHTML}`;
  box.querySelectorAll('.mpt').forEach(el => {
    el.style.cursor = 'pointer';
    el.addEventListener('click', () => location.href = `gene.html?id=${encodeURIComponent(el.dataset.sym)}`);
  });
}

function flagsHTML(flags, cfgFlags) {
  return '<div class="flagset">' + cfgFlags.map(f => {
    const on = flags[f.key];
    return `<span class="flag ${f.key} ${on ? 'on' : ''}" title="${f.desc}">${f.label}</span>`;
  }).join('') + '</div>';
}

function render() {
  const d = STATE.data;
  const sc = d.scores, fl = d.flags;
  // header
  let th = `<tr>
    <th data-key="symbol" class="${cls('symbol')}">Gene<span class="arrow">${arrow('symbol')}</span></th>
    <th data-key="protein_name" class="${cls('protein_name')}">Protein<span class="arrow">${arrow('protein_name')}</span></th>`;
  sc.forEach(s => th += `<th data-key="score:${s.key}" class="${cls('score:' + s.key)}" title="${s.desc}">${s.label}<span class="arrow">${arrow('score:' + s.key)}</span></th>`);
  th += `<th>Flags</th><th data-key="modality" class="${cls('modality')}">Modality<span class="arrow">${arrow('modality')}</span></th></tr>`;
  document.getElementById('thead').innerHTML = th;

  // rows
  let rows = d.genes.filter(g => {
    if (STATE.q) {
      const hay = (g.symbol + ' ' + g.protein_name).toLowerCase();
      if (!hay.includes(STATE.q)) return false;
    }
    for (const fk of STATE.flags) if (!g.flags[fk]) return false;
    return true;
  });
  rows.sort(cmp);
  const tb = rows.map(g => {
    let r = `<tr onclick="location.href='gene.html?id=${encodeURIComponent(g.symbol)}'">
      <td><span class="sym">${g.symbol}</span></td>
      <td class="pname">${g.protein_name || ''}</td>`;
    sc.forEach(s => r += `<td class="scorecell">${barHTML(g.scores[s.key], s)}</td>`);
    r += `<td>${flagsHTML(g.flags, fl)}</td>`;
    r += `<td><span class="modality-tag">${g.modality || 'Undetermined'}</span></td></tr>`;
    return r;
  }).join('');
  document.getElementById('tbody').innerHTML = tb || '<tr><td class="loading">No matches.</td></tr>';
  document.getElementById('count').textContent = `${rows.length} / ${d.genes.length} genes`;

  document.querySelectorAll('thead th[data-key]').forEach(th => th.onclick = () => setSort(th.dataset.key));
}

function val(g, key) {
  if (key.startsWith('score:')) { const v = g.scores[key.slice(6)]; return v === null ? -Infinity : v; }
  return (g[key] || '').toString().toLowerCase();
}
function cmp(a, b) {
  const va = val(a, STATE.sortKey), vb = val(b, STATE.sortKey);
  if (va < vb) return -1 * STATE.sortDir;
  if (va > vb) return 1 * STATE.sortDir;
  return 0;
}
function setSort(k) {
  if (STATE.sortKey === k) STATE.sortDir *= -1;
  else { STATE.sortKey = k; STATE.sortDir = k.startsWith('score:') ? -1 : 1; }
  render();
}
const cls = k => STATE.sortKey === k ? 'sorted' : '';
const arrow = k => STATE.sortKey === k ? (STATE.sortDir > 0 ? '▲' : '▼') : '↕';

function buildChips() {
  const box = document.getElementById('flagchips');
  box.innerHTML = STATE.data.flags.map(f =>
    `<span class="chip" data-flag="${f.key}" title="${f.desc}">${f.label}</span>`).join('');
  box.querySelectorAll('.chip').forEach(c => c.onclick = () => {
    const k = c.dataset.flag;
    if (STATE.flags.has(k)) { STATE.flags.delete(k); c.classList.remove('active'); }
    else { STATE.flags.add(k); c.classList.add('active'); }
    render();
  });
}

fetch('data/index.json').then(r => r.json()).then(d => {
  STATE.data = d;
  document.getElementById('site-title').textContent = d.site.title;
  document.getElementById('site-subtitle').textContent = d.site.subtitle;
  document.getElementById('pill-count').textContent = `${d.n_genes} genes`;
  if (window.renderFooter) window.renderFooter(d.site);
  document.getElementById('pill-note').textContent = '30-cytokine module';
  document.getElementById('search').addEventListener('input', e => { STATE.q = e.target.value.toLowerCase().trim(); render(); });
  buildChips();
  render();
  renderMatrix(d);
  window.addEventListener('themechange', () => renderMatrix(d));
}).catch(e => {
  document.getElementById('tbody').innerHTML = `<tr><td class="loading">Failed to load data.json — serve over HTTP (see README).</td></tr>`;
});
