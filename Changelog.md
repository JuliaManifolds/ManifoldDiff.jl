# Changelog

All notable changes to this Julia package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.4.3] June 24, 2025

### Added

* Compatibility bounds for `FiniteDifferences.jl` (0.12), `FiniteDiff.jl` (2.0), `ForwardDiff.jl` (0.10), `ReverseDiff.jl` (1.0), `Zygote.jl` (0.7), `DifferentiationInterface.jl` (0.7).

## [0.4.2] February 5, 2025

### Changed

* Increased ManifoldsBase.jl compatibility to include 1.0

## [0.4.1] December 2, 2024

### Added

* Basic usage example in docs.

## [0.4.0] November 27, 2024

### Changed

* Switch from manual backend creation to [ADTypes.jl](https://github.com/SciML/ADTypes.jl) + [DifferentiationInterface.jl](https://github.com/JuliaDiff/DifferentiationInterface.jl)
* Julia compat lower bound bumped from 1.6 to 1.10 (the new LTS)

## [0.3.13] November 13, 2024

### Added

* `jacobian_exp_argument` and its mutating variant, `jacobian_exp_argument!`.
* `jacobian_exp_basepoint` and its mutating variant, `jacobian_exp_basepoint!`.
* `jacobian_log_argument` and its mutating variant, `jacobian_log_argument!`.
* `jacobian_log_basepoint` and its mutating variant, `jacobian_log_basepoint!`.

## [0.3.12] September 5, 2024

### Added

* an individual logo that still resembles the `Manifolds.jl` family but also features a ∂.

## [0.3.11] August 28, 2024

### Changed

* Support for `Manifolds.jl` 0.10.

## [0.3.10] December 13, 2023

### Added

* Compatibility with `RecursiveArrayTools` v3.
* CI testing on Julia 1.10-RC.

## [0.3.9] - December 8, 2023

### Added

* proximal map of the distance function, `prox_distance`, and its mutating variant, `prox_distance!`
* this changelog
* A GitHub Action to check that the Changelog is updated on every PR.
