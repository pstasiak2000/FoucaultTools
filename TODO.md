# FoucaultTools TODO

This document tracks planned functionality and development tasks for the
FoucaultTools ecosystem.

---

# FoucaultBase

Core data structures and operations shared by the other FoucaultTools
packages.

## Fields

- [x] Implement basic `Field` operations
  - [x] Addition
  - [x] Subtraction
  - [x] Scalar multiplication
  - [x] Element-wise operations
  - [x] Other basic arithmetic operations

- [x] Implement basic vector-field operations
  - [x] Vector-field addition
  - [x] Vector-field subtraction
  - [x] Scalar multiplication
  - [x] Dot product
  - [ ] Cross product
  - [ ] Norm / magnitude

- [ ] Implement differential field operations
  - [ ] Gradient
  - [ ] Divergence
  - [ ] Curl
  - [ ] Laplacian

## Wavenumbers and Fourier transforms

- [ ] Develop `WaveNumbers` data structures for 1D, 2D and 3D
  - [x] Dimension-independent representation
  - [x] Complex FFT support
  - [x] Real FFT support
  - [ ] Add tests

- [ ] Add FFT planning structures
  - [ ] Forward FFT plans
  - [ ] Inverse FFT plans
  - [ ] Real-to-complex plans
  - [ ] Complex-to-real plans
  - [ ] Investigate plan reuse for repeated transforms

- [ ] Implement Fourier transforms for fields
  - [ ] Forward transform
  - [ ] Inverse transform
  - [ ] In-place transforms
  - [ ] Real FFT transforms

- [ ] Implement Fourier transforms for vector fields
  - [ ] Forward transform
  - [ ] Inverse transform
  - [ ] In-place transforms
  - [ ] Real FFT transforms

- [ ] Investigate efficient interaction between `Field`, `VectorField`
      and `WaveNumbers`

## I/O

- [ ] Implement field I/O
  - [ ] Read fields from files
  - [ ] Write fields to files
  - [ ] Support existing Foucault data formats
  - [ ] Document file formats
  - [ ] Add support for selecting data type when reading

- [ ] Investigate binary field formats for efficient I/O

---

# VFTools

Tools for working with vortex-filament data.

## Vortex data structures

- [ ] Design vortex data structures
  - [ ] Vortex line representation
  - [ ] Vortex filament representation
  - [ ] Vortex configuration / tangle representation
  - [ ] Periodic boundary representation
  - [ ] Vortex metadata

- [ ] Investigate efficient representations for large vortex tangles

## Vortex I/O

- [ ] Implement vortex data reading
  - [ ] Read individual vortex configurations
  - [ ] Read vortex tangles
  - [ ] Support existing Foucault vortex-file formats

- [ ] Implement vortex data writing
  - [ ] Write individual configurations
  - [ ] Write vortex tangles

## Basic vortex operations

- [ ] Vortex length
- [ ] Total vortex line length
- [ ] Vortex-line density
- [ ] Vortex tangent vectors
- [ ] Vortex curvature
- [ ] Vortex separation
- [ ] Vortex reconnection operations
- [ ] Periodic boundary handling
- [ ] Vortex statistics

---

# NSTools

Tools for analysing the normal-fluid / Navier–Stokes fields.

## Field quantities

- [ ] Calculate kinetic energy
- [ ] Calculate enstrophy
- [ ] Calculate energy spectra
- [ ] Calculate velocity statistics
- [ ] Calculate vorticity statistics
- [ ] Calculate other relevant field quantities

## Differential operators

- [ ] Derivatives
- [ ] Gradient
- [ ] Divergence
- [ ] Curl
- [ ] Laplacian
- [ ] Higher-order derivatives

- [ ] Implement derivatives in physical space
- [ ] Implement derivatives in Fourier space
- [ ] Compare physical-space and Fourier-space implementations

## Spectral analysis

- [ ] Fourier transforms of velocity fields
- [ ] Energy spectrum
- [ ] Enstrophy spectrum
- [ ] Wavenumber shell averaging
- [ ] Spectral derivatives
- [ ] Spectral filtering

## Statistics

- [ ] PDFs of field quantities
- [ ] Structure functions
- [ ] Velocity increment statistics
- [ ] Vorticity statistics
- [ ] Intermittency measures
- [ ] Higher-order moments

---

# General / Infrastructure

## Testing

- [ ] Add unit tests to `FoucaultBase`
- [ ] Add unit tests to `VFTools`
- [ ] Add unit tests to `NSTools`
- [ ] Add tests for 1D, 2D and 3D functionality
- [ ] Add tests for real and complex FFTs
- [ ] Add tests for periodic boundary conditions
- [ ] Add performance benchmarks for critical operations

## Documentation

- [ ] Document all public types
- [ ] Document all public functions
- [ ] Add examples for common workflows
- [ ] Set up `Documenter.jl`
- [ ] Create API reference
- [ ] Create user guide
- [ ] Add examples for typical FoucaultTools workflows

## Performance

- [ ] Identify allocation-heavy operations
- [ ] Add in-place (`!`) versions of expensive operations
- [ ] Benchmark FFT operations
- [ ] Investigate FFTW plan reuse
- [ ] Benchmark field operations
- [ ] Benchmark vortex operations
- [ ] Investigate multithreading where appropriate

## Package architecture

- [ ] Keep `FoucaultBase` independent of higher-level analysis packages
- [ ] Minimise dependencies between `VFTools` and `NSTools`
- [ ] Ensure common functionality lives in `FoucaultBase`
- [ ] Define clear public APIs for each package
- [ ] Review package interfaces as functionality grows
