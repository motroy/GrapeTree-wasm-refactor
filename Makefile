# Makefile for GrapeTree WebAssembly compilation
# Requires Emscripten SDK

CXX = em++
CXXFLAGS = -std=c++17 -O3 \
           -s WASM=1 \
           -s ALLOW_MEMORY_GROWTH=1 \
           -s MODULARIZE=1 \
           -s EXPORT_NAME='GrapeTreeWASMModule' \
           -s EXPORTED_RUNTIME_METHODS='["ccall","cwrap"]' \
           --bind \
           -I./src/cpp \
           -I./third_party/json/include

# Add single-file for easier distribution (optional)
SINGLE_FILE_FLAG = # -s SINGLE_FILE=1

SOURCES = src/cpp/wasm_interface.cpp

OUTPUT_DIR = build
OUTPUT_JS = $(OUTPUT_DIR)/grapetree.js
OUTPUT_WASM = $(OUTPUT_DIR)/grapetree.wasm

.PHONY: all clean test install-deps

all: $(OUTPUT_JS)

$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

$(OUTPUT_JS): $(SOURCES) | $(OUTPUT_DIR)
	$(CXX) $(CXXFLAGS) $(SINGLE_FILE_FLAG) $(SOURCES) -o $(OUTPUT_JS)
	@echo "✓ Build complete: $(OUTPUT_JS) and $(OUTPUT_WASM)"
	@ls -lh $(OUTPUT_DIR)

clean:
	rm -rf $(OUTPUT_DIR)
	@echo "✓ Build directory cleaned"

# Build optimized production version
production: CXXFLAGS += -Os --closure 1
production: clean $(OUTPUT_JS)
	@echo "✓ Production build complete (optimized)"

# Build debug version with source maps
debug: CXXFLAGS += -g -s ASSERTIONS=1 -s SAFE_HEAP=1 \
                   --source-map-base http://localhost:8080/build/
debug: clean $(OUTPUT_JS)
	@echo "✓ Debug build complete (with source maps)"

# Install nlohmann/json dependency
install-deps:
	@echo "Installing dependencies..."
	@mkdir -p third_party/json
	@if [ ! -d "third_party/json/include" ]; then \
		git clone --depth 1 https://github.com/nlohmann/json.git third_party/json-src && \
		cp -r third_party/json-src/include third_party/json/ && \
		rm -rf third_party/json-src && \
		echo "✓ nlohmann/json installed"; \
	else \
		echo "✓ nlohmann/json already installed"; \
	fi

# Run tests (requires Node.js)
test: $(OUTPUT_JS)
	@echo "Running tests..."
	node tests/test_suite.js

# Serve locally for testing
serve:
	@echo "Starting local server on http://localhost:8080"
	@cd web && python3 -m http.server 8080

# Build everything needed for deployment
deploy: production
	@echo "Preparing deployment..."
	@cp $(OUTPUT_JS) web/
	@cp $(OUTPUT_WASM) web/
	@echo "✓ Files copied to web/ directory"
	@echo "Deploy the contents of the web/ directory"

# Print build info
info:
	@echo "GrapeTree WebAssembly Build System"
	@echo "=================================="
	@echo "Compiler: $(CXX)"
	@echo "Flags: $(CXXFLAGS)"
	@echo "Output: $(OUTPUT_JS)"
	@echo ""
	@echo "Targets:"
	@echo "  make              - Build development version"
	@echo "  make production   - Build optimized version"
	@echo "  make debug        - Build with debugging"
	@echo "  make clean        - Remove build files"
	@echo "  make test         - Run test suite"
	@echo "  make serve        - Start local server"
	@echo "  make deploy       - Build and prepare for deployment"
	@echo "  make install-deps - Install dependencies"

# Check if Emscripten is available
check-emscripten:
	@which em++ > /dev/null || (echo "Error: Emscripten not found. Install from https://emscripten.org" && exit 1)
	@echo "✓ Emscripten found: $$(em++ --version | head -n1)"
