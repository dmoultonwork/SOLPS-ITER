# Up-down mirror symmetrisation tools for double-null grids

Python 3 (standard library only) tools used to make a B2.5 grid, its Eirene mesh and the Eirene surface/source inputs
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
`symmetrise_neutral_restart.py` hard-codes the masses of the nine virtual neutral species of the STEP case
(D2, D, ... and the He/Ar atoms and molecules, in that order) and the proton mass is the CODATA-2018 value (option
`--proton-mass` in `symmetrise_state.py`; it must match the mass used by the linked solver). These are STEP-specific and
should be generalised before the tools are used for other species sets.
