#!/bin/bash
# setup_onnxruntime.sh — Download ONNX Runtime for GNN phasing correction
#
# Run this before 'make' to enable the GNN module.
# Not needed if you don't use 'longphase gnn'.

set -e

VERSION="1.17.1"
ARCH=$(uname -m)

case "$ARCH" in
    x86_64)  PLATFORM="linux-x64" ;;
    aarch64) PLATFORM="linux-aarch64" ;;
    *)       echo "Unsupported architecture: $ARCH"; exit 1 ;;
esac

TARBALL="onnxruntime-${PLATFORM}-${VERSION}.tgz"
URL="https://github.com/microsoft/onnxruntime/releases/download/v${VERSION}/${TARBALL}"
DIR="onnxruntime-${PLATFORM}-${VERSION}"

if [ -d "onnxruntime" ] && [ -f "onnxruntime/lib/libonnxruntime.so" ]; then
    echo "ONNX Runtime already set up at onnxruntime/"
    echo "  $(ls onnxruntime/lib/libonnxruntime.so*)"
    exit 0
fi

echo "Downloading ONNX Runtime v${VERSION} for ${PLATFORM}..."
if command -v wget &>/dev/null; then
    wget -q --show-progress "$URL"
elif command -v curl &>/dev/null; then
    curl -LO "$URL"
else
    echo "ERROR: wget or curl required"
    exit 1
fi

echo "Extracting..."
tar xzf "$TARBALL"
ln -sf "$DIR" onnxruntime
rm -f "$TARBALL"

echo ""
echo "ONNX Runtime v${VERSION} ready at onnxruntime/"
echo "Now run: make -j \$(nproc)"
