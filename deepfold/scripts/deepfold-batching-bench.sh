#!/bin/bash


OUTPUT_FILE="outputs/bench/batchingbench.log"
PROOF_SIZE_FILE="batch.csv"

echo "Deepfold batching benchmarking..."

cargo bench -p batch -- --nocapture --quiet > $OUTPUT_FILE


if [ $? -eq 0 ]; then
    echo "Benchmark results have been written to $OUTPUT_FILE"
else
    echo "Benchmark failed to run"
fi

echo "Deepfold batch proofsizing..."
cargo test -p batch --release -- --nocapture --quiet
cp batch/batch.csv outputs/deepfold/$PROOF_SIZE_FILE
echo "Deepfold proofsize has been written to $PROOF_SIZE_FILE"



