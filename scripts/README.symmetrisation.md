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

Units: B2.5 geometry is in metres, Eirene `fort.33` in centimetres; check units before comparing.
The tools do not copy neutral restart data and do not symmetrise plasma states.
