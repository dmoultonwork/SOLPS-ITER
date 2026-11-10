# Up-down mirror symmetrisation tools for double-null grids

Python 3 tools (standard library only, except `symmetrise_neutral_restart.py` and `symmetrise_attached_side.py`, which also need
numpy and scipy) used to make a B2.5 grid, its Eirene mesh and the Eirene surface and source inputs exactly mirror symmetric, so
that any asymmetry seen later comes from the plasma evolution and not from the grid. Developed for the STEP double-null geometry;
not yet tested on other grids. Every tool states in its `--help` or docstring what it does **not** certify.

The Eirene input is the usual formatted input file (`fort.1`).

Typical order:

1. `audit_b2_reflection.py b2fgmtry LO:HI [LO:HI ...]`: read-only check of how far an existing `b2fgmtry` is from exact
   mirror symmetry for the given logical reflection blocks (inclusive guard-cell endpoints in x).
2. `symmetrise_carre.py --block LO:HI source.geo destination.geo`: write a reflected CARRE `.geo` input. Then regenerate
   `b2fgmtry` with the native `b2ag`, with the line `'b2agbs_average_flux' '1'` in `b2ag.dat` (option added in B2.5 by
   the same branch): it averages the reference flux over all valid western faces of a flux tube so that the metrics are
   mirror invariant.
3. `mirror_eirene_mesh.py --original DIR --old-b2 ... --new-b2 ... --output DIR`: build a symmetric candidate Eirene
   mesh from the verified lower half (moves plasma vertices with the B2 quadrilateral map, clips vacuum triangles,
   reflects, rebuilds neighbour links). It writes `fort.33`-`fort.35` and the reflection map of the cells
   (`triangle_reflection.txt`). Audit the result with
   `audit_mesh_symmetry.py --b2-geometry ... --eirene-vertices ... --eirene-triangles ... --output ...` (read-only, reports
   candidate reflected pairs and connectivity; nearest-point pairing is a candidate, never proof).
4. `symmetrise_eirene_inputs.py fort.1 pairing.txt new_fort.1`: mirror the end points of the additional surfaces and the
   coordinates of the surface sources in the formatted Eirene input file.
5. `test_mesh_audit.py`, `test_symmetrise_eirene_inputs.py`: small fixtures for the audit and the input tool.

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
The tools do not copy neutral restart data (except step 7) and do not symmetrise plasma states (except step 6).

Assumptions to check before use on another case: `audit_boundary_symmetry.py` expects 24 species per boundary entry,
`symmetrise_eirene_inputs.py` swaps the surface labels -4 <-> -3 and 5 <-> 6 by default (option `--swap A:B`, repeatable, to
change them), `symmetrise_neutral_restart.py` hard-codes the masses (2, 2, 2, 4, 4, 4, 4, 4, 4 proton masses) of the nine
virtual neutral species of the STEP restart file, in file order, and the proton mass is the CODATA-2018 value (option
`--proton-mass` in `symmetrise_state.py`; it must match the mass used by the linked solver). These are STEP-specific and
should be generalised before the tools are used for other species sets.

## Worked example

File names are illustrative. `old/` is the original case (`b2fgmtry`, `fort.33`/`fort.34`/`fort.35`, the formatted Eirene input
`fort.1`, a `b2fstate`, a `fort.13` and `b2.boundary.parameters`); `new/` receives the symmetric case. `LO1:HI1` and `LO2:HI2`
are the inclusive guard-cell endpoints in x of the two blocks that are mirror images of each other (together they must cover
-1..nx).

1. How far is the original grid from mirror symmetry?

       audit_b2_reflection.py old/b2fgmtry LO1:HI1 LO2:HI2

2. Reflect the CARRE geometry, then regenerate the B2 geometry with the native `b2ag`. The `b2ag.dat` in `new/` must contain
   the line `'b2agbs_average_flux' '1'` (this branch of B2.5). Repeat the audit on the new `b2fgmtry`:

       symmetrise_carre.py --block LO1:HI1 --block LO2:HI2 old/case.geo new/case.geo
       (run b2ag in new/)
       audit_b2_reflection.py new/b2fgmtry LO1:HI1 LO2:HI2

   On the STEP case this took the corner mismatch from 3e-4 m to 0, the cell-volume mismatch from 4.6% to 3e-14 and the
   `hx` mismatch from 6.8% to 1e-13.

3. Build the mirrored Eirene mesh (`new/mesh/fort.33`-`fort.35` and `new/mesh/triangle_reflection.txt`) and audit it:

       mirror_eirene_mesh.py --original old --old-b2 old/b2fgmtry --new-b2 new/b2fgmtry --output new/mesh
       audit_mesh_symmetry.py --b2-geometry new/b2fgmtry --eirene-vertices new/mesh/fort.33 \
           --eirene-triangles new/mesh/fort.34 --output new/mesh_audit.report

4. Mirror the additional surfaces and surface sources of the formatted Eirene input. `pairing.txt` is prepared by hand (these
   tools do not produce it): one line `i j` for each additional surface of block 3b (numbered 1..N in the file), giving its mirror
   partner, and `i i` for a surface that is its own mirror image; lines starting with `#` are comments. The mapping must be
   a complete involution:

       symmetrise_eirene_inputs.py old/fort.1 pairing.txt new/fort.1

   For every pair, the surface in the upper half receives exactly the mirror image of the end points of its partner, and the same
   is done for the source coordinates of the strata with a surface source; all other parameters of paired surfaces and sources
   must agree, or the script stops. The changed lines are listed. On the original STEP input this changed one additional-surface
   end point and the matching gas-puff source line by 0.001 cm, the only departures from mirror symmetry in the file.

5. Project the restart state and the neutral background onto the symmetric grid (the neutral tool reads
   `new/mesh/triangle_reflection.txt`):

       symmetrise_state.py --source old/b2fstate --old-b2 old/b2fgmtry --new-b2 new/b2fgmtry --output new/b2fstati
       symmetrise_neutral_restart.py old/fort.13 new/mesh new/fort.13

   On the STEP case the species inventories were conserved to 1e-16 and the electron and ion energies exactly.

6. Check the result and the boundary-condition input (each tool writes a report to the file named last):

       audit_state_symmetry.py new/b2fstati new/b2fgmtry new/state_audit.report
       measure_plasma_asymmetry.py new/b2fstati new/b2fgmtry new/asymmetry.report
       audit_boundary_symmetry.py old/b2.boundary.parameters new/boundary_audit.report

7. Run the case with the new `b2fgmtry` and `fort.30` (from `b2ag`), the mirrored `fort.33`-`fort.35`, the mirrored `fort.1`,
   `b2fstati` and `fort.13`. The surface and source physics are not certified by these tools: audit them separately.
