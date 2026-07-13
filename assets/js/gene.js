// Per-gene detail page: scores, Mol* structure, external links, Perturb-seq plots.
const PALETTE = {
  blue:['#5a3a2e','#d97757'], green:['#4a4230','#bfa98a'],
  amber:['#5a3f20','#e0a458'], red:['#4a2118','#cc6d52'],
};
const qs = new URLSearchParams(location.search);
const ID = qs.get('id');

function scoreBar(v, cfg) {
  if (v === null || v === undefined)
    return `<div class="barbox"><div class="bar"></div><span class="val na">n/a</span></div>`;
  const [c0, c1] = PALETTE[cfg.palette] || PALETTE.blue;
  const frac = Math.max(0, Math.min(1, (v - cfg.min) / (cfg.max - cfg.min)));
  return `<div class="barbox" title="${cfg.desc}">
    <div class="bar"><span style="width:${(frac*100).toFixed(1)}%;background:linear-gradient(90deg,${c0},${c1})"></span></div>
    <span class="val">${v.toFixed(2)}</span></div>`;
}

function effBars(items, title, note, axisMax) {
  const M = axisMax || 6;                       // global, fixed across all genes
  // axis scale ticks (−M .. +M) drawn once above the bars
  const ticks = [-M, -M/2, 0, M/2, M];
  const scale = ticks.map(t => {
    const pct = 50 + (t / M) * 50;              // 0 → 50%, +M → 100%, −M → 0%
    return `<span class="escl" style="left:${pct}%">${t>0?'+':''}${t}</span>`;
  }).join('');
  const rows = items.map(i => {
    const frac = Math.min(1, Math.abs(i.effect) / M);
    const up = i.effect >= 0;
    const w = (frac * 50).toFixed(1);
    const color = up ? 'var(--clay)' : 'var(--eff-down)';   // clay = induced by KD, slate = suppressed
    const style = up ? `left:50%;width:${w}%;background:${color}` : `right:50%;width:${w}%;background:${color}`;
    const sig = i.sig;
    return `<div class="eff${sig ? ' sig' : ''}" title="${i.name}: log₂FC ${i.effect>=0?'+':''}${i.effect.toFixed(2)} (Stim8hr)${sig ? ' — significant (FDR<0.05)' : ' — n.s.'}">
      <span class="nm">${i.name}${sig ? ' <b class="sigdot">●</b>' : ''}</span>
      <span class="track"><span class="mid"></span><span class="fill" style="${style};opacity:${sig?1:0.38}"></span></span>
      <span class="num"${sig?'':' style="opacity:.55"'}>${i.effect>=0?'+':''}${i.effect.toFixed(2)}</span></div>`;
  }).join('');
  return `<h2>${title}</h2>
    <div class="effscale"><span class="escl-lab" style="left:0">← suppressed by KD</span><span class="escl-lab" style="right:0;left:auto">induced by KD →</span></div>
    <div class="effscale effscale-ticks">${scale}</div>
    <div class="effbars">${rows}</div>
    <div class="subtle" style="margin-top:8px">${note}</div>`;
}

// Landscape: this gene highlighted among all 128 (x = selectivity, y = log2FC of representative cytokine)
function niceTicks(lo, hi, n) {
  const span = hi - lo || 1;
  const raw = span / n;
  const mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) || mag;
  const t0 = Math.ceil(lo / step) * step, out = [];
  for (let v = t0; v <= hi + 1e-9; v += step) out.push(Math.abs(v) < 1e-9 ? 0 : +v.toFixed(4));
  return out;
}

function landscapeSVG(points, gene) {
  if (!points || !points.length) return '';
  const W = 360, H = 320, padL = 52, padR = 16, padT = 14, padB = 46;
  const xs = points.map(p => p.sel), ys = points.map(p => p.lfc);
  const minX = 0, maxX = Math.max(...xs) * 1.03;
  const yPad = (Math.max(...ys) - Math.min(...ys)) * 0.04;
  const minY = Math.min(...ys) - yPad, maxY = Math.max(...ys) + yPad;
  const sx = x => padL + (x - minX) / (maxX - minX || 1) * (W - padL - padR);
  const sy = y => H - padB - (y - minY) / (maxY - minY || 1) * (H - padT - padB);
  const me = points.find(p => p.gene === gene);
  const dots = points.filter(p => p.gene !== gene).map(p =>
    `<circle cx="${sx(p.sel).toFixed(1)}" cy="${sy(p.lfc).toFixed(1)}" r="2.6"
      style="fill:var(--muted)" opacity="0.4"/>`).join('');
  // axis ticks with numeric values
  const xt = niceTicks(minX, maxX, 5), yt = niceTicks(minY, maxY, 5);
  const xticks = xt.map(t => `
    <line x1="${sx(t).toFixed(1)}" y1="${(H-padB).toFixed(1)}" x2="${sx(t).toFixed(1)}" y2="${(H-padB+4).toFixed(1)}" style="stroke:var(--muted)" stroke-width="0.8"/>
    <text x="${sx(t).toFixed(1)}" y="${(H-padB+15).toFixed(1)}" text-anchor="middle" class="ltick">${t}</text>`).join('');
  const yticks = yt.map(t => `
    <line x1="${(padL-4).toFixed(1)}" y1="${sy(t).toFixed(1)}" x2="${padL.toFixed(1)}" y2="${sy(t).toFixed(1)}" style="stroke:var(--muted)" stroke-width="0.8"/>
    <text x="${(padL-7).toFixed(1)}" y="${(sy(t)+3.3).toFixed(1)}" text-anchor="end" class="ltick">${t}</text>`).join('');
  const y0 = sy(0);
  let star = '', annot = '';
  if (me) {
    const cxp = sx(me.sel), cyp = sy(me.lfc);
    star = `<circle cx="${cxp.toFixed(1)}" cy="${cyp.toFixed(1)}" r="7"
      style="fill:var(--amber);stroke:var(--star-stroke)" stroke-width="1.6"><title>${gene}: selectivity ${me.sel.toFixed(2)}, log2FC ${me.lfc.toFixed(2)} on ${me.cyt}</title></circle>`;
    // callout label naming the representative cytokine, placed to avoid the top/right edges
    const lx = cxp > W * 0.6 ? cxp - 8 : cxp + 10;
    const anchor = cxp > W * 0.6 ? 'end' : 'start';
    const ly = cyp < padT + 24 ? cyp + 18 : cyp - 11;
    annot = `<text x="${lx.toFixed(1)}" y="${ly.toFixed(1)}" text-anchor="${anchor}" class="lcyt">${me.cyt}</text>`;
  }
  return `<h2>Position in the regulator landscape</h2>
    ${me ? `<div class="lcaption">Strongest effect on <b class="lcyt-inline">${me.cyt}</b> · <span class="ldir">${me.dir}</span> · log₂FC ${me.lfc>=0?'+':''}${me.lfc.toFixed(2)} · selectivity ${me.sel.toFixed(2)}</div>` : ''}
    <svg viewBox="0 0 ${W} ${H}" width="100%" style="max-height:320px">
      <line x1="${padL}" y1="${(H-padB).toFixed(1)}" x2="${(W-padR).toFixed(1)}" y2="${(H-padB).toFixed(1)}" style="stroke:var(--line)"/>
      <line x1="${padL}" y1="${padT}" x2="${padL}" y2="${(H-padB).toFixed(1)}" style="stroke:var(--line)"/>
      ${(minY < 0 && maxY > 0) ? `<line x1="${padL}" y1="${y0.toFixed(1)}" x2="${(W-padR).toFixed(1)}" y2="${y0.toFixed(1)}" style="stroke:var(--line)" stroke-dasharray="3 3"/>` : ''}
      ${xticks}${yticks}
      ${dots}${star}${annot}
      <text x="${((padL+(W-padR))/2).toFixed(1)}" y="${H-6}" text-anchor="middle" class="laxis">selectivity index →</text>
      <text x="13" y="${((padT+(H-padB))/2).toFixed(1)}" text-anchor="middle" transform="rotate(-90 13 ${((padT+(H-padB))/2).toFixed(1)})" class="laxis">log₂FC on ${me ? me.cyt : 'target cytokine'}</text>
    </svg>
    <div class="subtle">Gold marker = <b>${gene}</b> among all 128 regulators. x = selectivity (focused → right); y = log₂FC on its most strongly affected cytokine.</div>`;
}

// Pocket table beside the 3D viewer
function pocketTable(pockets) {
  if (!pockets || !pockets.length) return '';
  const rows = pockets.map(p =>
    `<tr><td>Pocket ${p.rank}</td>
      <td class="mono">${p.druggability!=null?p.druggability.toFixed(3):'—'}</td>
      <td class="mono">${p.volume!=null?Math.round(p.volume):'—'} Å³</td>
      <td class="mono">${p.residues.length}</td></tr>`).join('');
  return `<table class="pocktab">
    <thead><tr><th></th><th>Druggability</th><th>Volume</th><th>Residues</th></tr></thead>
    <tbody>${rows}</tbody></table>`;
}

function rnaDiagram(rna) {
  if (!rna || !rna.mrna_len) return '';
  const W = 640, H = 96, padL = 12, padR = 12, y = 44, barH = 16;
  const L = rna.mrna_len;
  const sx = p => padL + (p / L) * (W - padL - padR);
  const cds = (rna.cds_start && rna.cds_end);
  // backbone (full mRNA = UTR-colored), CDS block on top — colors via CSS vars (theme-aware)
  let svg = `<svg class="rna-diagram" viewBox="0 0 ${W} ${H}" width="100%" style="max-height:110px">`;
  svg += `<rect x="${sx(1).toFixed(1)}" y="${y}" width="${(sx(L)-sx(1)).toFixed(1)}" height="${barH}" rx="3" style="fill:var(--diag-utr)"/>`;
  if (cds) {
    svg += `<rect x="${sx(rna.cds_start).toFixed(1)}" y="${y-3}" width="${(sx(rna.cds_end)-sx(rna.cds_start)).toFixed(1)}" height="${barH+6}" rx="3" style="fill:var(--clay)"/>`;
    svg += `<text x="${((sx(rna.cds_start)+sx(rna.cds_end))/2).toFixed(1)}" y="${y+barH+18}" text-anchor="middle" style="fill:var(--accent)">CDS ${rna.cds_start}–${rna.cds_end}</text>`;
    svg += `<text x="${sx(1).toFixed(1)}" y="${y-9}" style="fill:var(--muted)">5′UTR</text>`;
    svg += `<text x="${sx(L).toFixed(1)}" y="${y-9}" text-anchor="end" style="fill:var(--muted)">3′UTR</text>`;
  }
  // scale ticks
  svg += `<text x="${sx(1).toFixed(1)}" y="${y+barH+18}">0</text>`;
  svg += `<text x="${sx(L).toFixed(1)}" y="${y+barH+18}" text-anchor="end">${L.toLocaleString()} nt</text>`;
  // real ASO / siRNA target sites — siRNA above the bar, ASO below, colored by modality
  (rna.sites || []).forEach(s => {
    const x = sx(s.pos);
    const isSi = s.modality === 'siRNA';
    const col = isSi ? 'var(--site)' : 'var(--eff-down)';
    const cy = isSi ? (y - 3) : (y + barH + 3);
    svg += `<circle class="site" cx="${x.toFixed(1)}" cy="${cy.toFixed(1)}" r="3.4"
      style="fill:${col};stroke:var(--site-ring)" stroke-width="1"><title>${s.modality} site @ ${s.start}–${s.end} nt · accessibility ${s.score!=null?s.score:'—'}\n${s.seq||''}</title></circle>`;
  });
  svg += `</svg>`;
  const nSi = (rna.sites||[]).filter(s=>s.modality==='siRNA').length;
  const nAso = (rna.sites||[]).filter(s=>s.modality==='ASO').length;
  return `<h2>mRNA map — ${rna.transcript || ''}</h2>${svg}
    <div class="rna-legend">
      <span class="k"><span class="sw" style="background:var(--diag-utr)"></span>UTR</span>
      <span class="k"><span class="sw" style="background:var(--clay)"></span>CDS</span>
      <span class="k"><span class="sw" style="background:var(--site);border-radius:50%;width:10px;height:10px"></span>siRNA site (${nSi})</span>
      <span class="k"><span class="sw" style="background:var(--eff-down);border-radius:50%;width:10px;height:10px"></span>ASO site (${nAso})</span>
    </div>
    <div class="subtle" style="margin-top:6px">Points mark qualifying nucleic-acid target windows (accessibility + efficacy + specificity gates). Hover for coordinates and sequence.</div>`;
}

function linkRow(links) {
  return `<h2>External resources</h2><div class="linkgrid">` +
    links.map(l => `<a href="${l.url}" target="_blank" rel="noopener">${l.label} ↗</a>`).join('') +
    `</div>`;
}

async function loadMolstar(g) {
  if (!g.structure) {
    document.getElementById('viewer').innerHTML =
      '<div class="loading">No AlphaFold model available.</div>';
    return;
  }
  try {
    const viewer = await molstar.Viewer.create('viewer', {
      layoutIsExpanded: false, layoutShowControls: false,
      layoutShowSequence: true, layoutShowLog: false,
      viewportShowExpand: true, viewportShowSelectionMode: false,
      pdbProvider: 'rcsb', emdbProvider: 'rcsb',
    });
    await viewer.loadStructureFromUrl(g.structure, 'pdb', false);
    // Color pocket-proximal residues via the occupancy channel (1.0 = residue lines a
    // druggable pocket, 0.0 = elsewhere). The occupancy theme's DEFAULT palette is a dark
    // "purples" scale — override it with a grey→clay ramp whose clay end is read from the
    // live --clay CSS var so it matches the caption swatch in BOTH themes, and re-apply on
    // theme change.
    const applyPocketTheme = async () => {
      try {
        const plugin = viewer.plugin;
        const structs = plugin.managers.structure.hierarchy.current.structures;
        const comps = structs.flatMap(s => s.components);
        const clay = cssHex('--clay', 0xD97757);      // live theme clay
        const grey = cssHex('--muted', 0xC9C2B6);     // live theme muted grey
        await plugin.managers.structure.component.updateRepresentationsTheme(comps, {
          color: 'occupancy',
          colorParams: { domain: [0, 1], list: { kind: 'interpolate', colors: [grey, clay] } },
        });
      } catch (themeErr) {
        console.warn('pocket theme not applied:', themeErr);  // structure still renders
      }
    };
    await applyPocketTheme();
    window.addEventListener('themechange', applyPocketTheme);
  } catch (e) {
    document.getElementById('viewer').innerHTML =
      '<div class="loading">Mol* failed to load structure (' + e + ').</div>';
  }
}

// Resolve a CSS custom property to a 0xRRGGBB int for Mol* (which wants numeric colors).
function cssHex(varName, fallback) {
  try {
    const v = getComputedStyle(document.documentElement).getPropertyValue(varName).trim();
    let m = v.match(/^#([0-9a-f]{6})$/i);
    if (m) return parseInt(m[1], 16);
    m = v.match(/^#([0-9a-f]{3})$/i);
    if (m) return parseInt(m[1].replace(/(.)/g, '$1$1'), 16);
    m = v.match(/rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)/i);
    if (m) return (+m[1] << 16) | (+m[2] << 8) | (+m[3]);
  } catch (e) {}
  return fallback;
}

fetch(`data/genes/${encodeURIComponent(ID)}.json`).then(r => {
  if (!r.ok) throw new Error('not found');
  return r.json();
}).then(g => {
  document.title = `${g.symbol} — Target Browser`;
  return Promise.all([
    fetch('data/index.json').then(r => r.json()),
    fetch('data/landscape.json').then(r => r.ok ? r.json() : {points:[]}).catch(() => ({points:[]})),
  ]).then(([idx, land]) => ({g, idx, land}));
}).then(({g, idx, land}) => {
  const sc = idx.scores;
  const dfields = idx.detail_fields || [];
  const cyto = (g.perturb && g.perturb.cytokines) || [];
  const detailHTML = dfields.length ? `<div class="detailgrid">` +
    dfields.map(d => `<div class="dcell"><span class="dk">${d.label}</span>
      <span class="dv">${g.details[d.key] || '—'}</span></div>`).join('') + `</div>` : '';
  const html = `
    <div class="gene-head">
      <span class="g">${g.symbol}</span>
      <span class="pn">${g.protein_name || ''}</span>
    </div>
    <div class="kv"><b>UniProt</b> ${g.uniprot||'—'} &nbsp;·&nbsp; <b>Length</b> ${g.length||'—'} aa
      &nbsp;·&nbsp; <b>Ensembl</b> ${g.ensembl_gene||'—'}
      &nbsp;·&nbsp; <b>Nominated modality</b> <span class="modality-tag">${g.modality||'Undetermined'}</span></div>
    ${g.function ? `<div class="func">${g.function}</div>` : ''}

    <div class="grid" style="margin-top:18px">
      <div class="panel">
        <h2>Druggability &amp; selectivity</h2>
        ${sc.map(s => `<div class="scorerow"><span class="lab" title="${s.desc}">${s.label}</span>
          <span class="barwrap">${scoreBar(g.scores[s.key], s)}</span></div>`).join('')}
        ${detailHTML}
        <div style="margin-top:16px">${linkRow(g.links)}</div>
      </div>
      <div class="panel">
        <h2>AlphaFold structure — druggable pockets</h2>
        <div id="viewer"></div>
        <div class="viewer-note">Drag to rotate · scroll to zoom. <b style="color:var(--clay)">Clay</b> = residues lining a fpocket-detected druggable pocket; grey = elsewhere. Model AF-${g.uniprot}-F1.</div>
        ${pocketTable(g.pockets)}
      </div>
    </div>

    ${g.rna ? `<div class="panel" style="margin-top:18px">${rnaDiagram(g.rna)}</div>` : ''}

    <div class="grid" style="margin-top:18px">
      <div class="panel">${landscapeSVG(land && land.points, g.symbol)}</div>
      <div class="panel">
        ${cyto.length
          ? effBars(cyto, 'Cytokine effects (all 30)',
              'Log₂ fold-change in Stim8hr on '+g.symbol+' knockdown, across all 30 cytokines in the module. Axis is fixed at ±'+((g.perturb&&g.perturb.axis_max)||6)+' for every gene so bar lengths are comparable target-to-target. <b class="sigdot">●</b> marks the '+((g.perturb&&g.perturb.n_sig)||cyto.filter(c=>c.sig).length)+' cytokine'+(((g.perturb&&g.perturb.n_sig)||0)===1?'':'s')+' this gene significantly regulates (FDR&lt;0.05); non-significant bars are dimmed.',
              (g.perturb && g.perturb.axis_max) || 6)
          : '<h2>Cytokine effects (all 30)</h2><div class="subtle">No cytokine data recorded for this gene.</div>'}
      </div>
    </div>
  `;
  document.getElementById('content').innerHTML = html;
  if (window.renderFooter) window.renderFooter(idx.site);
  loadMolstar(g);
}).catch(e => {
  document.getElementById('content').innerHTML =
    `<div class="loading">Gene "${ID}" not found. <a href="index.html">Back to list</a></div>`;
});
