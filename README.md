# Stochastic LDPM in Julia

A Julia 1.12 implementation of the Lattice Discrete Particle Model (LDPM) for fracture and compressive failure of ordinary and pervious concrete. It is a port and extension of Dr. Alessandro Fascetti's MATLAB code. The stochastic version (S-LDPM) adds a spatially correlated random field (SCRF) for the pore structure of pervious concrete.

The repository is built stage by stage. Every stage has fixed inputs and outputs and a set of checks it must pass before the next stage builds on it.

> **Status:** private, under construction (github.com/fascetti/juldpm). Stages 1–2 are in; the rest are plans. See [Status](#status).

---

## Units

All files use **mm, N, s, tonne, MPa**. Every file header repeats this.

| Quantity | Unit | Typical value |
|---|---|---|
| Length | mm | — |
| Stress, modulus | MPa (N/mm²) | E ≈ 45 000 MPa |
| Mass | tonne | — |
| Density | tonne/mm³ | 2.4e-9 for 2400 kg/m³ (**not** 2.4e-6, which is kg/mm³) |
| Fracture energy | N/mm | 0.04 N/mm = 40 N/m |
| Time step | s | 1e-5 |

Check every conversion; the density factor of 1000 has already caused one error in the legacy code.

---

## Pipeline

Seven stages, each in its own folder under `src/`:

| Stage | Folder | Owner | Produces |
|---|---|---|---|
| 1 Mesher | `src/mesher/` (+ `src/geometry/`) | Hasan (STL, box, cylinder); Radia, Baixi (FreeCAD) | `Nodes`, `Diameters` (0 = surface node), `Node_Face` (surface patch per node) |
| 2 Delaunay and clean-up | `src/delaunay/` | Hasan | `Tets` (Int32, TetGen); renumbered nodes, orphans removed |
| 3 Facets | `src/facets/` | Hasan | per-edge elements by default (D1): facet areas, centroids, local axes, lengths |
| 4 Compatibility, material, mass | `src/kinematics/` | Alessandro (B), Hasan | B operators, per-element material arrays, lumped mass, Δt |
| 5 Boundary conditions | `src/boundary/` | Radia | per-DOF tag and imposed velocity |
| 6 Solvers | `src/solvers/` | Alessandro, Hasan (explicit); Baixi (static, arc-length) | displacements, forces, reactions, element strains and stresses |
| 7 Output | `src/output/` | Hasan | force–displacement curve, VTK frames, timing report |

Porosity (`src/porosity/`) is optional and comes in later. The SCRF acts at stage 1, on particle placement, and at stages 3–4, through the effective-area factor ζ_e (Fascetti et al. 2022, Eq. 22).

The exact interface of every stage (array names, types, sizes, units) is in [`docs/STAGE_INTERFACES.md`](docs/STAGE_INTERFACES.md).

---

## Data layout: flat arrays only

The time loop works only on flat arrays, never on arrays of structs. Rewriting the multithreaded CPU solver this way cut its run time from 171 s to 5 s on the same machine (user-reported; mesh and step count to be recorded in `docs/BENCHMARKS.md`). The same layout runs on the GPU.

The rules:

- **One plain array per field.** Element `i` is index `i` in every array. No per-element struct, and no `Vector` stored inside a struct.
- **Fixed-size small vectors are `SVector`s,** for example the six B rows of an element. Indices are `Int32`; all real numbers are `Float64`.
- **Hot loops are functions that take every array as an argument.** No globals inside a loop, and no variable that is reassigned and captured by `@threads`.
- **No allocation inside the time loop.** Allocate once before it and update in place.
- **Assembly is a CSR gather.** Each DOF sums its own list of contributions, in ascending element order. That makes it race-free and deterministic, and bit-for-bit equal to a serial loop.
- **Node arrays are n × 3,** one row per node as in the MATLAB code. A hot loop may work on a transposed copy.

Pre-processing keeps the same layout: per-entity data are flat arrays (CSR for variable-length lists such as triangles per grid cell), bundled in a NamedTuple where a function needs many of them (the surface, `STL`). There are no per-particle, per-tet or per-triangle structs anywhere.

---

## Repository layout

```text
.
├── Project.toml, Manifest.toml   # the Julia environment (not a package)
├── includes.jl                   # ordered include list, one function per file
├── src/
│   ├── core/          # input reader, RNG, timing, gate printing
│   ├── geometry/      # surface: STL reading, box/cylinder generation, queries, cell list
│   ├── mesher/        # stage 1
│   ├── delaunay/      # stage 2
│   ├── facets/        # stage 3
│   ├── kinematics/    # stage 4
│   ├── boundary/      # stage 5
│   ├── solvers/       # stage 6: explicit/, static/
│   ├── porosity/      # SCRF, ζ_e (planned)
│   └── output/        # stage 7
├── examples/          # run_<test>_<specimen>.jl drivers
├── inputs/            # mesh_*.txt input files; geometry/ holds the STL files
├── test/              # runtests.jl; regression references later
├── install_packages.jl
└── docs/              # STAGE_INTERFACES.md, DECISIONS.md, NAMES.md, BENCHMARKS.md, PORT_NOTES.md
```

The code keeps the MATLAB-style include idiom: one function per file, included in order by `includes.jl`, with drivers in `examples/`. It is a Julia environment, not a package. Packaging can come later if the code goes to JOSS or the General registry.

---

## Getting started

```bash
git clone https://github.com/fascetti/juldpm.git
cd juldpm
julia install_packages.jl                      # once: installs, checks TetGen, writes Manifest.toml
julia --project=. -t 12 examples/run_mesh.jl inputs/mesh_box_75x75x150.txt
julia --project=. -t 4 test/runtests.jl
```

- `-t 12` sets the number of threads; use your core count.
- Every run prints each gate, the wall time of each stage, and the total. The first run in a session includes compilation.
- `seed = 0` in the input file draws a new seed and prints it; the seed is written to `mesh_nodes.csv`, so the same seed reproduces the same specimen.
- Inputs for a cylinder, the dogbone, a notched beam and a tube (hole test) are in `inputs/`.

---

## Conventions

- **Comment tags.** Every intentional difference from the MATLAB or legacy Julia code is tagged:
  - `FAITHFUL-NOTE`: kept as in the original, with the reason;
  - `CHANGED`: changed on purpose, with the reason;
  - `DEPRECATED`: kept only for comparison;
  - `FAITHFUL-BUG-x`: a known bug kept deliberately, numbered.
- **Timing.** Every solver prints its total wall time, plus a per-phase breakdown.
- **Randomness.** Pass `rng::AbstractRNG` explicitly and write the seed to the output. Never rely on the global RNG.
- **File writes.** Write output to a configured folder, never into `src/`.
- **Names.** Canonical variable names are in `docs/NAMES.md`. Rename a variable in one commit that updates every file using it.
- **Legacy code.** The uniaxial and three-point-bending legacy trees are frozen, read-only references, kept outside this repository. Nothing is copied from them; their two `types.jl` files diverged and must not be merged by copying.
- **Physics first.** A change that alters results states the size of the change and its physical cause in the pull request.

---

## Checks (gates)

Each stage ships with checks that run before its output is used, and print pass or fail:

| Stage | Checks |
|---|---|
| 1 Mesher | surface valid (watertight, outward normals, mm); particle–wall and particle–particle clearance; volume fraction against target; surface nodes on the surface; same seed gives identical nodes |
| 2 Delaunay | total tet volume equals specimen volume; no flat or outside tets; no orphan nodes |
| 3–4 | rigid-body motion gives zero facet strain; uniform strain gives the exact normal strain; assembled mass without scaling equals ρV |
| 6 Explicit | threaded equals serial bit for bit; quasi-static checks: R_top + R_bot ≈ 0 and E_kin ≪ W_ext before the peak |
| Flatten | CSR gather equals serial scatter bit for bit; flat element response equals the reference for every element |

The external benchmark is Lale et al. (2026): shared particle, tetrahedron and facet files, public on Zenodo (record 18197080), which can be read in place of our own mesher.

---

## Status

| Stage | State |
|---|---|
| 1 Mesher (STL, box, cylinder) | in the repo with gates; logic checked in a Python mirror; first Julia run pending |
| 1 Mesher (FreeCAD) | in progress; writes `mesh_nodes.csv` in the same format |
| 2 Delaunay (TetGen) | in the repo with gates; first Julia run pending |
| 3 Facets | next: per edge (D1), tet-sweep connectivity |
| 4 Compatibility, material, mass | planned |
| 5 Boundary conditions | planned |
| 6 Explicit solver | flat CPU prototype exists, one element per Delaunay edge; per-facet version waits on D1 |
| 6 Static solver | planned |
| 7 Output | planned |
| Porosity | planned, after stage 4 |
| GPU backend | planned, after the CPU solver is fixed |

---

## Decisions

All decisions, settled and open, are in [`docs/DECISIONS.md`](docs/DECISIONS.md). Settled so far: per-edge elements by default with a per-facet option (D1), TetGen for every geometry (D4), a required `pore_check` when porosity is on (D5), a named wall rule (D9), seeded StableRNG (D10), and sampled surface nodes instead of the legacy grid (D11). Still open: ε_V (D2), the surface-node spacing to use for comparison with Paper 1 (D11), loaded area on curved specimens (D13), and licence and release (D8).

---

## Contributing

- `main` is protected. Work on a short-lived branch, for example `mesher/stl-face-tags`, and open a pull request.
- Every pull request needs one approval and cannot be approved by its author.
- Each pull request states:
  - which stage it touches;
  - units checked;
  - the gates it ran and their output;
  - whether results change, by how much and why;
  - any use of AI tools.
- Keep pull requests small; one stage or one fix each.

---

## References

- Cusatis, G., Pelessone, D., Mencarelli, A. (2011). Lattice discrete particle model (LDPM) for failure behavior of concrete. I: Theory. *Cement and Concrete Composites* 33(9), 881–890.
- Fascetti, A., Bolander, J. E., Nisticò, N. (2018). Lattice discrete particle modeling of concrete under compressive loading: multiscale experimental approach for parameter determination. *Journal of Engineering Mechanics* 144(8), 04018058.
- Fascetti, A., et al. (2022). Stochastic lattice discrete particle modeling of fracture in pervious concrete. *Computer-Aided Civil and Infrastructure Engineering* 37(14).
- Lale, E., et al. (2026). arXiv:2603.13190. Benchmark data: Zenodo record 18197080.

---

## People

| Name | Role |
|---|---|
| Dr. Alessandro Fascetti | principal investigator; compatibility operators, explicit solver |
| Md Khairul Hasan | STL mesher, Delaunay, facets, explicit solver, output |
| Radia | FreeCAD mesher, boundary conditions |
| Baixi | FreeCAD mesher, static (arc-length) solver |

## Licence

Not yet chosen. The repository stays private until the UCF disclosure and the Pitt release are settled.
