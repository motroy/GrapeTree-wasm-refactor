#!/usr/bin/env node
// Render a GrapeTree image from the command line by driving the web app
// (GitHub Pages by default) in a headless browser. Requires: npm i playwright
import { chromium } from 'playwright';
import { parseArgs } from 'node:util';
import { writeFileSync } from 'node:fs';
import { resolve, extname } from 'node:path';

const { values: o } = parseArgs({
  options: {
    profile: { type: 'string' }, metadata: { type: 'string' }, hiercc: { type: 'string' },
    method: { type: 'string', default: 'MSTreeV2' },
    'hc-level': { type: 'string' }, 'shade-clusters': { type: 'boolean', default: false },
    'color-by': { type: 'string' }, scale: { type: 'string', default: '2' },
    out: { type: 'string', default: 'tree.svg' },
    url: { type: 'string', default: 'https://motroy.github.io/GrapeTree-wasm-refactor/' },
    example: { type: 'boolean', default: false }, help: { type: 'boolean', short: 'h' },
  },
});
if (o.help || !(o.profile || o.example)) {
  console.log(`Usage: grapetree-render.mjs --profile FILE [--metadata FILE] [--hiercc FILE]
         [--method MSTreeV2|MSTree] [--hc-level HC10] [--shade-clusters]
         [--color-by COLUMN|HierCC] [--out tree.svg|tree.png] [--scale 2]
         [--url BASE_URL] [--example]`);
  process.exit(o.help ? 0 : 1);
}
const ext = extname(o.out).toLowerCase();
if (ext !== '.svg' && ext !== '.png') throw new Error('--out must end in .svg or .png');

const browser = await chromium.launch();
try {
  const page = await browser.newPage({ viewport: { width: 1400, height: 1000 } });
  page.on('pageerror', e => console.error('page error:', e.message));
  await page.goto(o.url);
  await page.waitForFunction(() => typeof d3 !== 'undefined');

  if (o.example) {
    await page.click('#load-example-btn');
  } else {
    await page.setInputFiles('#file-input', resolve(o.profile));
    if (o.metadata) await page.setInputFiles('#metadata-input', resolve(o.metadata));
  }
  if (o.hiercc) await page.setInputFiles('#hiercc-input', resolve(o.hiercc));
  await page.waitForTimeout(1000);

  await page.selectOption('#method-select', o.method);
  await page.click('#compute-btn');
  await page.waitForSelector('#tree-canvas svg');

  if (o['hc-level'] || o['shade-clusters'] || o['color-by'] === 'HierCC') {
    await page.click('button[data-tab=clusters]');
    if (o['hc-level']) await page.selectOption('#hc-level-select', o['hc-level']);
    if (o['shade-clusters']) await page.check('#shade-clusters-check');
  }
  if (o['color-by'] && o['color-by'] !== 'HierCC') {
    await page.click('button[data-tab=style]');
    await page.selectOption('#color-select', o['color-by']);
  }

  // Wait for the force layout to settle, then fit the SVG to the tree's bounding box
  let prev = '';
  for (let stable = 0; stable < 3;) {
    await page.waitForTimeout(500);
    const b = await page.evaluate(() => {
      const r = document.querySelector('#tree-canvas svg > g').getBBox();
      return [r.x, r.y, r.width, r.height].map(Math.round).join(',');
    });
    stable = b === prev ? stable + 1 : 0;
    prev = b;
  }
  const svg = await page.evaluate(() => {
    const s = document.querySelector('#tree-canvas svg').cloneNode(true);
    const g = s.querySelector(':scope > g');
    const r = document.querySelector('#tree-canvas svg > g').getBBox();
    g.removeAttribute('transform');
    const pad = 30, w = Math.ceil(r.width + 2 * pad), h = Math.ceil(r.height + 2 * pad);
    s.setAttribute('viewBox', `${r.x - pad} ${r.y - pad} ${w} ${h}`);
    s.setAttribute('width', w); s.setAttribute('height', h);
    s.setAttribute('style', 'background: #fff');
    return new XMLSerializer().serializeToString(s);
  });

  if (ext === '.svg') {
    writeFileSync(o.out, svg);
  } else {
    const p2 = await (await browser.newContext({ deviceScaleFactor: Number(o.scale) })).newPage();
    await p2.setContent(`<body style="margin:0">${svg}</body>`);
    await (await p2.$('svg')).screenshot({ path: o.out });
  }
  console.log(`Wrote ${o.out}`);
} finally {
  await browser.close();
}
