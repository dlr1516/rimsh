# RIMSH — Rotation of Real Spherical Harmonics via the Lie algebra so(3)

MATLAB implementation accompanying our paper, developed at
RIMLab (Robotics and Intelligent Machines Laboratory), University of Parma.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24556000.svg)](https://doi.org/10.5281/zenodo.24556000)

> **If you use our work, please cite this DOI: <https://doi.org/10.5281/zenodo.24556000>**

The library computes the rotation matrices of **real spherical harmonics (RSH)** expansions from the
representation of the Lie algebra so(3) in the RSH basis. For a rotation of angle θ about the unit
axis **u** (Rodrigues vector **x** = θ**u**), the matrix of band *l* is

$$
M_l(\mathbf{x}) = \exp\big(x_1 L_x^{(l)} + x_2 L_y^{(l)} + x_3 L_z^{(l)}\big) = F_l(\mathbf{u})\, M_{l,z}(\theta)\, F_l(\mathbf{u})^\top ,
$$

where $L_x^{(l)}, L_y^{(l)}, L_z^{(l)}$ are sparse real generators, $F_l(\mathbf{u})$ is a real orthogonal
matrix that depends only on the axis, and $M_{l,z}(\theta)$ is the elementary rotation about *z*
(entries $\cos m\theta$ on the diagonal, $\pm\sin m\theta$ on the anti-diagonal).
This **anti-diagonal factorization** (called *Lie-AD* in the paper) avoids matrix exponentials and
complex eigen-decompositions, gives matrices that are orthogonal to round-off, and makes rotations
about a fixed axis, closed-form derivatives and Newton alignment on SO(3) inexpensive.

## Requirements

- MATLAB R2016b or later (no toolbox needed), or GNU Octave ≥ 6 for the library and the tests.
- The paper script `experiments/run_paper_experiments.m` uses local functions in a script and
  `timeit`, so it needs MATLAB.

## Quick start

```matlab
setup_path                                   % adds src/ to the path

lmax = 16;
x    = (pi/3) * [1; 1; 0] / sqrt(2);         % 60 deg about [1 1 0]/sqrt(2)
[Ml, M] = rsh_rot_lie_eig(lmax, x);          % M{l+1}: (2l+1)x(2l+1) matrix of band l

c     = von_mises_fisher_rsh_coeffs(lmax, [0;0;1], 30);   % RSH coefficients of a vMF density
c_rot = rsh_rot_coeffs(x, c);                % rotated expansion (x or a 3x3 matrix)
rsh_plot(c_rot, 90, 'rotated vMF');
```

A longer walk-through is in [`experiments/example_basic_usage.m`](experiments/example_basic_usage.m).

## Reproducing the paper

```matlab
run('experiments/run_paper_experiments.m')
```

The script runs every experiment of the paper (convention checks, timing and accuracy vs *l*,
rotation of von Mises–Fisher distributions, sweeps about a fixed axis, closed-form derivatives and
Newton alignment), writes the figures to `figures/` and the numbers to
`results/rsh_rot_experiments_results.mat`. Random draws are seeded (`rng(3)`); timings depend on the
machine. The full run (*l* up to 64) takes from several minutes to about an hour; reduce
`cfg.lList` for a quick check.

| Figure | File |
|---|---|
| Computation time vs *l* | `figures/rsh_time.pdf` |
| Discrepancy from Ivanic–Ruedenberg and orthogonality | `figures/rsh_accuracy.pdf` |
| Rotation of von Mises–Fisher distributions | `figures/rsh_vmf.pdf` |
| Many angles about a fixed axis | `figures/rsh_fixed_axis.pdf` |
| Closed-form derivatives and Newton alignment | `figures/rsh_derivatives.pdf` |

## Tests

```matlab
run('tests/run_tests.m')
```

Checks the commutation relations of the generators, agreement with `expm` and with the
Ivanic–Ruedenberg recursion, orthogonality, equivariance $Y_l(R\mathbf{v}) = M_l(R)\,Y_l(\mathbf{v})$,
the homomorphism property, rotation of closed-form vMF coefficients, gradient/Hessian against
finite differences and Newton alignment. Runs in a few seconds in MATLAB and Octave.

## Repository layout

```
src/                         core library
  rsh_lie_generator.m        generators Lx, Ly, Lz of band l
  rsh_rot_lie_eig.m          rotation matrices of bands 0..lmax, anti-diagonal factorization (proposed)
  rsh_rot_to_z.m             orthogonal factor F_l(u) of the factorization
  rsh_rot_lie.m              same matrices via expm (reference implementation)
  rsh_rot_coeffs.m           rotate a packed RSH coefficient vector
  rsh_eval.m, csh_eval.m     evaluate real / complex spherical harmonics
  sh_lm_to_index.m           packed index of (l, m)
  von_mises_fisher_rsh_coeffs.m   closed-form RSH coefficients of a vMF density
  rsh_plot.m                 lat/lon plot of an RSH expansion
  baselines/                 Ivanic–Ruedenberg recursion (rsh_rot_ivanic and helpers)
  registration/              correlation, gradient/Hessian, initial guess and Newton
                             alignment on SO(3) (experimental)
  utils/                     exp/log maps of SO(3) and small helpers
experiments/
  run_paper_experiments.m    reproduces all experiments and figures of the paper
  example_basic_usage.m      short tutorial
  demo_von_mises_fisher.m    vMF function vs its truncated RSH series
tests/run_tests.m            self-checks
figures/                     figures of the paper (generated)
```

## Conventions

- Coefficients of bands 0..*l*<sub>max</sub> are packed band after band, *m* = −*l*..*l*, at index
  *l*² + *l* + *m* + 1 (`sh_lm_to_index`); total length (*l*<sub>max</sub>+1)².
- Real SH are orthonormal on the sphere and built from MATLAB's Schmidt semi-normalized
  `legendre(l, cos θ, 'sch')`, with sin(*m*φ) for *m* < 0 and cos(*m*φ) for *m* > 0.
- Rotations are passed as Rodrigues vectors **x** = θ**u** or 3×3 matrices *R* = exp([**x**]<sub>×</sub>);
  `rsh_rot_ivanic(R, l)` and `rsh_rot_lie_eig(l, log_SO3(R))` return the same matrices.
- Known limitation: for axes at a tiny but non-zero angle ρ from ±*z*, the accuracy of `rsh_eval`
  limits *F<sub>l</sub>* to about eps/ρ (see the note in `rsh_rot_to_z.m`).

## Citation

If you use our work, please cite this DOI: **<https://doi.org/10.5281/zenodo.24556000>**
(see also [`CITATION.cff`](CITATION.cff); GitHub shows a "Cite this repository" button from it).

```bibtex
@software{rimlab_rimsh_2026,
  author    = {{RIMLab (Robotics and Intelligent Machines Laboratory), University of Parma}},
  title     = {{RIMSH: Rotation of Real Spherical Harmonics via the Lie algebra so(3)}},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24556000},
  url       = {https://doi.org/10.5281/zenodo.24556000}
}
```

## References

- J. Ivanic, K. Ruedenberg, "Rotation matrices for real spherical harmonics. Direct determination by
  recursion", *J. Phys. Chem.* 100(15), 6342–6347, 1996; Additions and corrections, *J. Phys. Chem. A*
  102(45), 9099–9100, 1998.

## License

Released under the [MIT License](LICENSE).
