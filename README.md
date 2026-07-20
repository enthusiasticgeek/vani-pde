# vani-pde

Finite-difference PDE solver library for the [vāṇी compiler](https://github.com/enthusiasticgeek/vani-compiler).

Depends on [vani-matrix](https://github.com/enthusiasticgeek/vani-matrix) for
`mat_zeros`/`mat_solve`, used by the two elliptic (steady-state) solvers.
Does **not** depend on vani-calculus (see "Design decisions" below).

## Add to your project

```toml
# vani.toml
[deps]
pde = { registry = "kosh", version = "^0.1" }
```

```sh
vanic add pde
vanic build
```

## What's included (v0.1.0 — complete; see TODO.md)

| Equation | Dimension | Functions |
|---|---|---|
| Laplace/Poisson (elliptic) | 1D | `pde_laplace1d_assemble`, `pde_laplace1d_rhs`, `pde_laplace1d_solve` |
| Laplace/Poisson (elliptic) | 2D | `pde_laplace2d_assemble`, `pde_laplace2d_rhs`, `pde_laplace2d_solve` |
| Heat (parabolic) | 1D | `pde_heat1d_stability_r`, `pde_heat1d_step`, `pde_heat1d_solve` |
| Heat (parabolic) | 2D | `pde_heat2d_stability_r`, `pde_heat2d_step`, `pde_heat2d_solve` |
| Wave (hyperbolic) | 1D | `pde_wave1d_courant`, `pde_wave1d_first_step`, `pde_wave1d_step`, `pde_wave1d_solve` |
| Wave (hyperbolic) | 2D | `pde_wave2d_courant`, `pde_wave2d_first_step`, `pde_wave2d_step`, `pde_wave2d_solve` |
| Grid utility | -- | `pde_grid2d_index` |

## Design decisions (read this before using the library)

This package has the widest design surface of anything in the Kosh math
ecosystem so far, so the choices below were made explicitly rather than
left implicit in the code:

- **Grid encoding**: flat row-major `Vec<f64>`, explicit dimension args, same
  convention as vani-matrix/vani-tensor. A 1D grid is `Vec<f64>` of length
  `n`. A 2D grid is `Vec<f64>` of length `nx*ny`; point `(i, j)` is at flat
  offset `i*ny + j` (`pde_grid2d_index`) -- the same row-major layout as a
  vani-matrix matrix or a vani-tensor rank-2 tensor.
- **Boundary conditions**: Dirichlet only. 1D solvers take two scalars
  (`left_bc`, `right_bc`). 2D solvers take a `boundary: Vec<f64>` the same
  length as the grid, whose border entries are used and whose interior
  entries are ignored (avoids a second, differently-indexed array).
  Neumann/periodic BCs are not implemented -- see TODO.md.
- **Elliptic equations** (Laplace/Poisson, time-independent) are solved
  *directly*: the finite-difference stencil is assembled into a dense
  matrix (boundary rows become identity rows) and solved once via
  vani-matrix's `mat_solve`. Dense assembly costs `O((nx*ny)^2)` memory --
  fine at this library's "modest grid size" scope, the same tradeoff
  vani-geometry made choosing `O(n^2)` convex hull over a divide-and-conquer
  algorithm.
- **Parabolic (heat) and hyperbolic (wave) equations** are solved via
  **explicit** time marching (FTCS for heat, central-difference for wave) --
  the simplest correct scheme, chosen over implicit/Crank-Nicolson methods
  to keep v0.1.0 tractable. Explicit schemes are only *conditionally*
  stable. This library computes the stability number for you
  (`pde_heat1d_stability_r`, `pde_heat2d_stability_r`, `pde_wave1d_courant`,
  `pde_wave2d_courant`) but does **not** enforce it at runtime -- check it
  yourself before calling `_step`/`_solve` (heat needs the result `<= 0.5`,
  wave needs it `<= 1.0`), same caller-trust convention as the rest of this
  ecosystem.
- **No vani-calculus dependency**, despite the roadmap sketch listing it:
  no natural reuse point was found without forcing an awkward fit.
  vani-calculus's ODE solvers step a scalar state via `fn(f64) -> f64`;
  PDE time-marching advances a whole-grid `Vec<f64>` state per step, a
  different shape of problem. This package depends on vani-matrix only.

## Correctness

Every solver is validated against a closed-form solution, not just
plausible-looking output: 1D/2D Laplace against a harmonic polynomial
(`x^2 - y^2`, exact to floating-point precision since the 5-point
Laplacian is exact for quadratics), 1D/2D heat against
`sin(pi x) exp(-alpha pi^2 t)`-style solutions, and 1D/2D wave against
`sin(pi x) cos(pi c t)`-style solutions. See `tests/`.

## What this library does NOT provide

These are already vāṇी compiler builtins — call them directly, no import needed:

`sin` `cos` `exp` `sqrt` `abs` `f64_pi()` `push` `pop` `len` `set` `vec`

## License

MIT
