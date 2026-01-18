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

## Documentation

See the full [refactoring plan](grapetree-wasm-refactor-plan.md) for architecture details and implementation roadmap.

## License

GNU GPL v3.0 - Same as original GrapeTree
