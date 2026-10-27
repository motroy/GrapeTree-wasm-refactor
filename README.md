# GrapeTree WebAssembly Edition

A high-performance WebAssembly port of [GrapeTree](https://github.com/achtman-lab/GrapeTree), enabling phylogenetic tree visualization and analysis directly in the browser with near-native performance.

## Overview

This is a complete refactoring of GrapeTree's core algorithms from Python to C++ and WebAssembly, providing:

- **5-10x faster** tree generation compared to the Python implementation
- **Zero installation** - runs entirely in modern web browsers
- **Offline capable** - can work without internet connection
- **Cross-platform** - same performance on all operating systems
- **Scalable** - handles datasets with 100,000+ strains efficiently

## Quick Start

### 1. Prerequisites

Install [Emscripten](https://emscripten.org/docs/getting_started/downloads.html):

```bash
git clone https://github.com/emscripten-core/emsdk.git
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh
```

### 2. Build

```bash
cd build
make release
```

### 3. Run

```bash
python3 -m http.server 8000
# Open http://localhost:8000/src/web/index.html
```

## Features

- ✅ MSTree & MSTreeV2 algorithms
- ✅ Distance matrix calculations
- ✅ Newick format export
- ✅ Interactive web interface
- 🚧 Neighbor-Joining (in progress)

## Command-line rendering

`tools/grapetree-render.mjs` drives the web app (the GitHub Pages deployment by default) in a headless browser and saves the tree as SVG or PNG:

```bash
cd tools && npm install && npx playwright install chromium
node grapetree-render.mjs --profile data.profile --metadata meta.txt \
  --color-by Country --out tree.svg
node grapetree-render.mjs --profile data.profile --hc-level HC10 \
  --shade-clusters --color-by HierCC --out tree.png --scale 3
```

Options: `--method`, `--hiercc FILE`, `--hc-level`, `--shade-clusters`, `--color-by`, `--scale`, `--url` (e.g. a local server), `--example`. The legend is not included in the output.

## Documentation

See the full [refactoring plan](grapetree-wasm-refactor-plan.md) for architecture details and implementation roadmap.

## License

GNU GPL v3.0 - Same as original GrapeTree
