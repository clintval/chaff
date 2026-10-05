# chaff

[![Build Status](https://github.com/clintval/chaff/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/clintval/chaff/actions/workflows/ci.yml?query=branch%3Amain)
[![Coverage Status](https://coveralls.io/repos/github/clintval/chaff/badge.svg?branch=main)](https://coveralls.io/github/clintval/chaff?branch=main)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-rust-dea588.svg)](https://www.rust-lang.org/)

Separate somatic variant calls from library-preparation damage artifacts.

## Introduction

Duplex and UMI sequencing suppress most sequencing and amplification errors, but not damage that library preparation copies onto both strands of a molecule before the strands are tagged.
The tool `chaff` scores each somatic call by where its alternate allele sits on the molecules that carry it, compared with the molecules that carry the reference allele at the same site.

Install from source:

```bash
cargo install --git https://github.com/clintval/chaff
```

## Development and Testing

See the [contributing guide](./CONTRIBUTING.md) for more information.

> [!NOTE]
> Claude Code was used substantially in the development of `chaff`, most notably for ideation support, prototyping, and code generation.
> Until a v1 release, treat this project as AI-enabled and under active review.
