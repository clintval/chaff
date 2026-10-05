# Contributing

## Toolchain

Use [`pixi`](https://pixi.sh/) for a reproducible environment including a Rust toolchain.

## Building

```bash
cargo build --release
./target/release/chaff --help
```

The `dist` profile (`cargo build --profile dist`) creates a release optimized, stripped binary for distribution.

## Mandatory Pre-submission Checks

All must pass before a change is submitted:

```bash
cargo build --release
cargo test --all --no-fail-fast --verbose
cargo fmt -- --check
cargo clippy --all-targets --all-features -- -D warnings
```

Equivalent pixi tasks:

```bash
pixi run build
pixi run test
pixi run fmt-check
pixi run lint
```
