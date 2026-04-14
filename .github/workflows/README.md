                    Finite Element Discretization Library
                                   __
                       _ __ ___   / _|  ___  _ __ ___
                      | '_ ` _ \ | |_  / _ \| '_ ` _ \
                      | | | | | ||  _||  __/| | | | | |
                      |_| |_| |_||_|   \___||_| |_| |_|

                               https://mfem.org

This directory contains the GitHub CI scripts for MFEM.

The current CI workflows are:

## `test.yml`

Runs the Zig 0.16.0 build on a small matrix of serial library configurations.

- Debug static build
- ReleaseSafe static build
- ReleaseSafe shared build
- ReleaseSafe single-precision build

## `repo-check.yml`

Runs a number of static repository-level sanity checks.

- `file-headers-check` checks copyright date, license, etc. using the `--copyright`, `--license` and `--release` options of the `config/githooks/pre-push` script.

- `code-style` checks the code style using the `--style` option of the `config/githooks/pre-push` script.

- `documentation` checks the documentation build using the `tests/scripts/documentation` script.

- `branch-history` guards against accidental commits of large files using the `--history` option of the `config/githooks/pre-push` script.

## `codeql-analysis.yml`

Runs CodeQL analysis after a manual Zig build.

## `stale.yml`

Marks stale issues and pull requests according to the project policy.
