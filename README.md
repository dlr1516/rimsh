# RIMSH — RIMLab Spherical Harmonics
#### Copyright (C) 2026 Dario Lodi Rizzini, Ernesto Fontana, Francesco Pasetti


OVERVIEW
-------------------------------------------------

MATLAB utilities for manipulating **Spherical Harmonics** and, in particular, for
**rotating Real Spherical Harmonics (RSH)**.

A real spherical harmonic series of degree $L$ is a vector of $(L+1)^2$ coefficients
expanded on the real basis $Y_{lm}(\theta,\varphi)$. Rotating such a series is the
operation that links the representation of a function to the frame it is expressed in:
given a rotation $R \in SO(3)$, the rotated coefficients are $\mathbf{c}' = \mathbf{R}(R)\,\mathbf{c}$.
RIMSH provides the evaluation of the basis, the construction of $\mathbf{R}(R)$, and the
tools needed to use RSH series inside optimisation and estimation problems
(gradients, Hessians, initial guesses, and Newton steps on $SO(3)$).

The rotation matrices follow a **Lie-algebraic formulation**: the generators of the
$so(3)$ representation acting on the degree-$L$ RSH subspace are derived from the
differential operators of infinitesimal rotations on the sphere, and $\mathbf{R}(R)$
is obtained through the exponential map. As an alternative that avoids a
general-purpose matrix exponential, $\mathbf{R}(R)$ is factorised into anti-diagonal
form — the rotation angle enters only through an elementary rotation about the
$z$ axis, while an orthogonal similarity matrix depends only on the rotation axis.
The formulation is free of the singularities affecting Euler-angle methods, remains
accurate to round-off at high degrees, and makes rotation-dependent functions of RSH
series explicitly differentiable. An implementation of the classical
Ivanic–Ruedenberg recursion is included for comparison.

### Reference

If you use this library in academic work, please cite:

> E. Fontana, D. Lodi Rizzini, F. Pasetti,
> *A Lie-Algebraic Formulation of Real Spherical Harmonic Rotation*,
> Zenodo, 2026. DOI: [10.5281/zenodo.24556000](https://zenodo.org/records/24556000)
> — [PDF](https://zenodo.org/records/24556000/files/2026_rsh_rot_arxiv.pdf?download=1)

```bibtex
@misc{rimsh2026,
  title     = {A Lie-Algebraic Formulation of Real Spherical Harmonic Rotation},
  author    = {Fontana, Ernesto and Lodi Rizzini, Dario and Pasetti, Francesco},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24556000},
  url       = {[doi.org](https://doi.org/10.5281/zenodo.24556000})
}
```


MAIN FUNCTIONS
-------------------------------------------------

### General mathematical utilities

Functions shared by the real (RSH) and complex (CSH) cases, plus the
linear-algebra helpers used throughout.

| Function | Description |
| --- | --- |
| `sh_lm_to_index` | Map a degree/order pair $(l,m)$ to its index in the coefficient vector (and back). |
| `sh_rot_index` | Index bookkeeping for the blocks of the RSH rotation matrix. |
| `legendre_poly` | Evaluate the Legendre polynomial $P_l(x)$ by recurrence. |
| `tchebychev_to_legendre_coeffs` | Convert Chebyshev coefficients into Legendre coefficients. |
| `pnebi_eval_0` | Evaluate the associated Legendre functions $P_l^0$ used to build the RSH basis. |
| `pnebi_eval_n` | Evaluate the associated Legendre functions $P_l^m$ of general order. |
| `exp_SO3` | Exponential map: rotation vector in $so(3)$ → rotation matrix in $SO(3)$. |
| `log_SO3` | Logarithm map: rotation matrix in $SO(3)$ → rotation vector in $so(3)$. |
| `adjoint_rodrigues` | Adjoint action of $SO(3)$ on $so(3)$ in Rodrigues form. |
| `factor_up` | Rank-one factorisation helper used by the anti-diagonal decomposition. |
| `get_opt` | Read an option from an options struct with a default value (shared by the optimisation routines). |
| `resetRands` | Reset the random number generator seed for reproducible experiments. |

### RSH evaluation and rotation

| Function | Description |
| --- | --- |
| `rsh_eval` | Evaluate the real spherical harmonic basis $Y_{lm}(\theta,\varphi)$ up to degree $L$. |
| `rsh_plot` | Plot an RSH series on the sphere. |
| `rsh_lie_generator` | Generator of the $so(3)$ representation acting on the degree-$L$ RSH subspace. |
| `rsh_rot_coeffs` | Rotate a set of RSH coefficients by a rotation matrix. |
| `rsh_rot_lie` | RSH rotation matrix from the Lie-algebraic formulation (exponential map). |
| `rsh_rot_lie_eig` | Similarity matrix from the eigenvectors of the RSH rotation matrix, by a re-orthogonalised ladder recurrence. |
| `rsh_rot_to_z` | RSH rotation matrix aligned with the $z$ axis (the similarity block of the factorisation). |
| `rsh_rot_to_z5` | Same as `rsh_rot_to_z`, closed-form case up to degree 5. |
| `rsh_rot_ivanic` | RSH rotation matrix by the classical Ivanic–Ruedenberg recursion (reference implementation). |
| `rsh_rot_ivanic_coeffs` | Rotate RSH coefficients with the Ivanic–Ruedenberg recursion. |
| `csh_eval` | Evaluate the complex spherical harmonic basis. |
| `csh_rot_wigner` | Rotation of complex spherical harmonics through Wigner $D$-matrices. |

### Use of RSH series: optimisation and modelling

Higher-level routines built on the rotation machinery.

| Function | Description |
| --- | --- |
| `rsh_corr_grad_hess` | Correlation between RSH series together with its closed-form gradient and Hessian on $SO(3)$ (via the Lie-algebraic linearisation). |
| `rsh_corr_init_guess` | Closed-form initial guess for the rotation aligning two RSH series (first-order approximation). |
| `rsh_corr_init_guess_2nd` | Second-order refinement of the initial guess. |
| `rsh_corr_ub_l` | Upper bound on the correlation for a given degree $L$, used to schedule or truncate the optimisation. |
| `rsh_opt_newton_so3` | Newton(-like) optimisation of a rotation on $SO(3)$ using `rsh_corr_grad_hess`. |
| `rsh_project_to_so3` | Project a matrix onto $SO(3)$ (nearest rotation). |
| `von_mises_fisher_rsh_coeffs` | RSH coefficients of a von Mises–Fisher distribution on the sphere (synthetic RSH series for tests and demos). |


DEPENDENCIES
-------------------------------------------------

