# DeepSPRINT

_Note: This is a prototype implementation and not ready for real applications!_

This repository contains a prototype implementation for the DeepSPRINT Isogeny Proof of Knowledge / Signature scheme. Running `cargo run --release` will run keygen, prove and verify for all parameter settings. First for a conjectured optimization and then for the prover variant. This will output the proof sizes.

Running `cargo bench` will run all benchmarks, both for inidvidual functions as well as the full keygen/prove/verify algorithms. Running `cargo bench --bench my_benchmark` will only run proving/verification and keygen benchmarks.

## DeepFold
The implementation of DeepFold comes from the official benchmarking repository. We have added prime fields to support our parameter sets.
