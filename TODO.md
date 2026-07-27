# vani-pde — TODO

> Compiler builtins that already exist and must NOT be reimplemented:
> `sin` `cos` `exp` `sqrt` `abs` `f64_pi()` `push` `pop` `len` `set` `vec`
>
> Depends on vani-matrix (`mat_zeros`, `mat_solve`) -- v0.1.0's only Kosh
> dependency. Does NOT depend on vani-calculus (see README's "Design
> decisions" for why, despite the roadmap sketch listing it).

---

## v0.1.0 — Implemented ✓

### Grid utility (1 function)
- [x] `pde_grid2d_index` -- row-major 2D flatten, matching
      vani-matrix/vani-tensor's convention

### 1D Laplace/Poisson, elliptic (3 functions)
- [x] `pde_laplace1d_assemble`, `pde_laplace1d_rhs`, `pde_laplace1d_solve`
      -- dense tridiagonal-ish assembly + vani-matrix's `mat_solve`.
      Validated against the exact solution of `u'' = -2` (`u(x) = x - x^2`)

### 2D Laplace/Poisson, elliptic (3 functions)
- [x] `pde_laplace2d_assemble`, `pde_laplace2d_rhs`, `pde_laplace2d_solve`
      -- 5-point stencil, dense `(nx*ny)^2` assembly + `mat_solve`.
      Validated against a harmonic polynomial (`x^2 - y^2`) matching the
      exact solution at every interior point, not just approximately --
      the 5-point Laplacian is exact for quadratics, a strong test

### 1D heat, parabolic (3 functions)
- [x] `pde_heat1d_stability_r`, `pde_heat1d_step`, `pde_heat1d_solve` --
      explicit FTCS. Validated against `sin(pi x) exp(-alpha pi^2 t)`

### 2D heat, parabolic (3 functions)
- [x] `pde_heat2d_stability_r`, `pde_heat2d_step`, `pde_heat2d_solve` --
      explicit FTCS. Validated against
      `sin(pi x) sin(pi y) exp(-alpha * 2 * pi^2 * t)`

### 1D wave, hyperbolic (4 functions)
- [x] `pde_wave1d_courant`, `pde_wave1d_first_step`, `pde_wave1d_step`,
      `pde_wave1d_solve` -- explicit central-difference, Taylor-series
      first step from `(u0, v0)`. Validated against `sin(pi x) cos(pi c t)`
      and the `steps=0` degenerate case (returns `u0` unchanged)

### 2D wave, hyperbolic (4 functions)
- [x] `pde_wave2d_courant`, `pde_wave2d_first_step`, `pde_wave2d_step`,
      `pde_wave2d_solve` -- same scheme in 2D. Validated against
      `sin(pi x) sin(pi y) cos(pi c sqrt(2) t)`

### Tests and examples
- [x] `tests/test_laplace.vani` -- 1D exact-polynomial check, 2D
      exact-harmonic-polynomial check (interior matches to 1e-9, not just
      "close")
- [x] `tests/test_heat.vani` -- 1D and 2D against closed-form decaying-sine
      solutions, plus explicit stability-number assertions
- [x] `tests/test_wave.vani` -- 1D and 2D against closed-form
      standing-wave solutions, the `steps=0` edge case, plus explicit
      Courant-number assertions
- [x] `examples/heat_diffusion_demo.vani` -- a hot spot on a 1D rod
      diffusing and smoothing over time
- [x] `examples/membrane_wave_demo.vani` -- a plucked 2D drum membrane,
      center-point displacement oscillating as the wave reflects off the
      clamped edges

### Safety annotations
- [x] `#[bounded_stack(bytes=N)]` on every function, budgets set to `vanic
      check`'s exact reported worst-case (largest: `pde_laplace2d_solve` at
      648 bytes, since it composes assemble + rhs + `mat_solve`'s own
      chain through `mat_zeros`)
- [x] No recursion anywhere in this library -- all time-marching is done
      with `while` loops calling `_step` repeatedly from `_solve`

---

## v0.1.4 (2026-07-27)

- [x] `pde_heat1d_step_neumann`/`pde_heat1d_solve_neumann` -- zero-flux
      (insulated) boundary conditions for the 1D heat equation, via the
      standard ghost-point trick (`du/dx=0` at an end means the ghost
      point mirrors its neighbor, e.g. `u[-1]=u[1]`, substituted into the
      same interior FTCS formula). Added as new functions alongside the
      existing Dirichlet `pde_heat1d_step`/`_solve`, not a BC-mode flag
      on them -- matches this package's existing one-function-per-
      PDE/BC-combination style (e.g. the separate `pde_laplace1d_*`/
      `pde_laplace2d_*` rather than a dimension flag).
- [x] `pde_heat1d_step_periodic`/`pde_heat1d_solve_periodic` -- periodic
      boundary conditions (the domain wraps: `u[-1]=u[n-1]`,
      `u[n]=u[0]`), so both ends use the OTHER end's real neighbor
      instead of a ghost/BC value.
- [x] `tests/test_heat.vani` extended with two closed-form checks: `du/dx=0`
      Neumann with `u(x,0)=cos(pi x)` (exact `cos(pi x) exp(-alpha pi^2 t)`,
      same eigenvalue magnitude as the existing Dirichlet sine case, same
      `1e-3` tolerance), and periodic with `u(x,0)=cos(2 pi x)` (exact
      `cos(2 pi x) exp(-4 alpha pi^2 t)`, a *looser* `5e-3` tolerance --
      confirmed empirically that this doubled-frequency mode genuinely
      has ~3x the FTCS truncation error at the same dx/dt, not just
      widened until the assertion happened to pass). Full suite +
      `vanic audit-safety` re-verified on both backends. No unrelated
      WCET/stack drift found in this package.
- Scoped to the 1D heat equation only in this pass (the simplest,
      clearest case) -- 2D heat and both 1D/2D wave equations still only
      support Dirichlet BCs; extending them the same way is a natural,
      contained follow-up, not attempted here to keep this change
      reviewable.

## Future

No v0.2.0 is currently planned. Candidates if a concrete need shows up:
Neumann/periodic BCs for the 2D heat equation and both wave equations
(1D done above; 2D heat and wave still Dirichlet-only), implicit or
Crank-Nicolson schemes for the heat equation (unconditionally stable,
would need vani-matrix's `mat_solve` per time step -- more expensive per
step but removes the `r <= 0.5` restriction), sparse matrix assembly for
the elliptic solvers instead of dense `O((nx*ny)^2)` (would need a new
`vani-sparse` package or an extension to vani-matrix, see kosh-index's
ROADMAP.md gap analysis), and higher-order (RK4-style) time stepping for
the wave equation.
