# Up-down mirror symmetrisation tools for double-null grids

Python 3 tools (standard library only, except `symmetrise_neutral_restart.py` and `symmetrise_attached_side.py`, which also need
numpy and scipy) used to make a B2.5 grid, its Eirene mesh and the Eirene surface/source inputs
exactly mirror symmetric, so that any asymmetry seen later comes from the plasma evolution and not from the grid.
Developed for the STEP double-null geometry; not yet tested on other grids. Every tool states in its `--help` or
docstring what it does **not** certify.

Typical order:

1. `audit_b2_reflection.py b2fgmtry LO:HI [LO:HI ...]`: read-only check of how far an existing `b2fgmtry` is from exact
   mirror symmetry for the given logical reflection blocks (inclusive guard-cell endpoints in x).
2. `symmetrise_carre.py --block LO:HI source.geo destination.geo`: write a reflected CARRE `.geo` input. Then regenerate
   `b2fgmtry` with the native `b2ag`, with the line `'b2agbs_average_flux' '1'` in `b2ag.dat` (option added in B2.5 by
   the same branch): it averages the reference flux over all valid western faces of a flux tube so that the metrics are
   mirror invariant.
3. `mirror_eirene_mesh.py --original fort.33... --old-b2 ... --new-b2 ... --output ...`: build a symmetric candidate Eirene
   mesh from the verified lower half (moves plasma vertices with the B2 quadrilateral map, clips vacuum triangles,
   reflects, rebuilds neighbour links). Audit the result with
   `audit_mesh_symmetry.py --b2-geometry ... --eirene-vertices ... --eirene-triangles ... --output ...` (read-only, reports
   candidate reflected pairs and connectivity; nearest-point pairing is a candidate, never proof).
4. `symmetrise_eirene_inputs.py reference pairing destination mesh`: mirror the surface and source coordinates of the
   Eirene JSON inputs, preserving surface IDs and physics; label permutations are reported separately.
5. `test_mesh_audit.py`: small analytic fixtures for the audit (`python3 test_mesh_audit.py`).

Restart (initial state) tools, used after the grid is symmetric:

6. `symmetrise_state.py`: project a restart onto mirror symmetry using conserved cell totals (species inventories, electron
   thermal energy, combined ion thermal plus parallel kinetic energy; the forbidden even part of the parallel momentum is
   removed; cached face fluxes are reflected with their face index shift). `symmetrise_attached_side.py` is the alternative
   start that copies one half onto the other instead of averaging the two halves.
7. `symmetrise_neutral_restart.py`: the same projection for the short virtual-neutral background of an Eirene restart
   (`fort.13`); rate caches are invalidated so Eirene rebuilds them.
8. `audit_state_symmetry.py`, `measure_plasma_asymmetry.py` and `audit_boundary_symmetry.py`: read-only checks of the
   reflected parities of a state, of the volume-weighted inventories and mirror residuals, and of the explicitly supplied
   boundary-condition values under the logical reflection.

Units: B2.5 geometry is in metres, Eirene `fort.33` in centimetres; check units before comparing.
The tools do not copy neutral restart data and do not symmetrise plasma states.

Assumptions to check before use on another case: `audit_boundary_symmetry.py` expects 24 species per boundary entry,
`symmetrise_eirene_inputs.py` hard-codes the surface label swaps -4 <-> -3 and 5 <-> 6,
`symmetrise_neutral_restart.py` hard-codes the masses (2, 2, 2, 4, 4, 4, 4, 4, 4 proton masses) of the nine virtual
neutral species of the STEP restart file, in file order, and the proton mass is the CODATA-2018 value (option
`--proton-mass` in `symmetrise_state.py`; it must match the mass used by the linked solver). These are STEP-specific and
should be generalised before the tools are used for other species sets.

## Worked example

File names are illustrative. `old/` is the original case (`b2fgmtry`, `fort.33`/`fort.34`/`fort.35`, the JSON Eirene inputs
`eirene.input.json`, `eirene_add_surfaces.input.json`, `eirene_physics_model.input.json`, a `b2fstate`, a `fort.13` and
`b2.boundary.parameters`); `new/` receives the symmetric case. `LO1:HI1` and `LO2:HI2` are the inclusive guard-cell endpoints in x
of the two blocks that are mirror images of each other (together they must cover -1..nx).

1. How far is the original grid from mirror symmetry?

       audit_b2_reflection.py old/b2fgmtry LO1:HI1 LO2:HI2

2. Reflect the CARRE geometry, then regenerate the B2 geometry with the native `b2ag`. The `b2ag.dat` in `new/` must contain
   the line `'b2agbs_average_flux' '1'` (this branch of B2.5). Repeat the audit on the new `b2fgmtry`:

       symmetrise_carre.py --block LO1:HI1 --block LO2:HI2 old/case.geo new/case.geo
       (run b2ag in new/)
       audit_b2_reflection.py new/b2fgmtry LO1:HI1 LO2:HI2

   On the STEP case this took the corner mismatch from 3e-4 m to 0, the cell-volume mismatch from 4.6% to 3e-14 and the
   `hx` mismatch from 6.8% to 1e-13.

3. Build the mirrored Eirene mesh (`fort.33`-`fort.35` and `mesh_generation.json`) and audit it (the audit lists nearest reflected
   points and the connectivity; it is a check, not proof of correct pairing):

       mirror_eirene_mesh.py --original old --old-b2 old/b2fgmtry --new-b2 new/b2fgmtry --output new/mesh
       audit_mesh_symmetry.py --b2-geometry new/b2fgmtry --eirene-vertices new/mesh/fort.33 \
           --eirene-triangles new/mesh/fort.34 --output new/mesh_audit.json

4. Mirror the surface and source coordinates of the JSON Eirene inputs. `pairing.json` is prepared by hand and is not produced by
   these tools: it holds `{"candidate_pairs": {"1": j1, "2": j2, ...}}`, the mirror partner of every additional surface
   (numbered 1..N), and must be a complete involution. `old` holds the three JSON files:

       symmetrise_eirene_inputs.py old pairing.json new/inputs new/mesh

5. Project the restart state and the neutral background onto the symmetric grid (the neutral tool reads
   `new/mesh/mesh_generation.json`):

       symmetrise_state.py --source old/b2fstate --old-b2 old/b2fgmtry --new-b2 new/b2fgmtry --output new/b2fstati
       symmetrise_neutral_restart.py old/fort.13 new/mesh new/fort.13

   On the STEP case the species inventories were conserved to 1e-16 and the electron and ion energies exactly.

6. Check the result and the boundary-condition input:

       audit_state_symmetry.py new/b2fstati new/b2fgmtry new/state_audit.json
       measure_plasma_asymmetry.py new/b2fstati new/b2fgmtry new/asymmetry.json
       audit_boundary_symmetry.py old/b2.boundary.parameters new/boundary_audit.json

7. Run the case with the new `b2fgmtry` and `fort.30` (from `b2ag`), the mirrored `fort.33`-`fort.35`, the mirrored JSON inputs,
   `b2fstati` and `fort.13`. The surface and source physics are not certified by these tools: audit them separately.
