# Up-down mirror symmetrisation tools for double-null grids

Python 3 tools (standard library only, except `symmetrise_neutral_restart.py` and `symmetrise_attached_side.py`, which also need
numpy and scipy) used to make a B2.5 grid, its Eirene mesh and the Eirene surface and source inputs exactly mirror symmetric, so
that any asymmetry seen later comes from the plasma evolution and not from the grid. Developed for the STEP double-null geometry;
not yet tested on other grids. Every tool states in its `--help` or docstring what it does **not** certify.

The Eirene input is the usual formatted input file (`fort.1`).

Typical order:

0. `find_reflection_blocks.py b2fgmtry`: read-only. Checks that the grid is a connected double null (two cuts at the same
   radial position, two poloidal blocks) and that each block is its own mirror image to within a tolerance (default a quarter
   of a cell size), and prints the block endpoints `LO:HI` that the following tools need. It exits with an error if the grid is
   not suitable, because there is then no point in symmetrising it. The classification follows the B2.5 cut conventions and
   was checked on STEP grids only.
1. `audit_b2_reflection.py b2fgmtry LO:HI [LO:HI ...]`: read-only check of how far an existing `b2fgmtry` is from exact
   mirror symmetry for the given logical reflection blocks (inclusive guard-cell endpoints in x).
2. `symmetrise_carre.py --block=LO:HI source.geo destination.geo`: write a reflected CARRE `.geo` input. Then regenerate
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
`symmetrise_eirene_inputs.py` infers which surface labels the mirror exchanges from the pairing and prints them (STEP: -4 <-> -3
and 5 <-> 6); the option `--swap A:B`, repeatable, overrides this. Check that the printed exchange is what you expect. `symmetrise_neutral_restart.py` hard-codes the masses (2, 2, 2, 4, 4, 4, 4, 4, 4 proton masses) of the nine
virtual neutral species of the STEP restart file, in file order, and the proton mass is the CODATA-2018 value (option
`--proton-mass` in `symmetrise_state.py`; it must match the mass used by the linked solver). These are STEP-specific and
should be generalised before the tools are used for other species sets.

## Worked example

All commands below are typed in a **tcsh** shell, one line each (do not break them with a backslash: a line break lost in copy and
paste silently drops the `--` of an option), from the `baserun` directory of your case. Set up once per shell session:

    source setup.csh            # in the top directory of the SOLPS-ITER tree: puts the tools of scripts/ on your PATH
    source SETUP/mpi            # b2run then uses the MPI executables; also 'source SETUP/openmp' if you built with OpenMP
    cd runs/<your case>/baserun
    mkdir old new               # old/ gets copies of your original files (below), new/ receives the symmetric case

File names are illustrative. `old/` is the original case: the CARRE grid file `case.geo` (the `.geo` file that `carre` wrote and
`b2ag` reads), the `b2ag.dat` used to build the grid, `b2fgmtry`, `fort.33`/`fort.34`/`fort.35`, the formatted Eirene input
`fort.1`, a `b2fstate`, a `fort.13` and `b2.boundary.parameters`. `new/` receives the symmetric case. `LO1:HI1` and `LO2:HI2`
are the inclusive guard-cell endpoints in x of the two mirror blocks (together they must cover -1..nx); step 0 finds them.
Copy the files listed for `old/` from your original `baserun` and run directory into `old/` (the tools only read `old/`). The starting
state is `b2fstate` if you take the end of a run, or `b2fstati` if you take an input state: use whichever you have, under its own name in
the commands below. Your original `baserun` also holds `b2ah.dat`, `b2ai.dat`, `b2ar.dat` and `b2frates`, which step 7 needs.

0. Is the grid a connected double null, and which are the blocks?

       find_reflection_blocks.py old/b2fgmtry

   It prints, for example for the STEP grid, `poloidal blocks: -1:72, 73:194`, the topology and the largest mismatch of
   mirrored cell centres, and suggests the commands of the next steps.

1. How far is the original grid from mirror symmetry?

       audit_b2_reflection.py old/b2fgmtry LO1:HI1 LO2:HI2

2. Reflect the CARRE geometry, then regenerate the B2 geometry with the native `b2ag`. For `b2ag.dat`, take the *same*
   `old/b2ag.dat` that built the original grid, append the line `'b2agbs_average_flux' '1'` (option of this branch of B2.5),
   and save it as `new/b2ag.dat`; change nothing else in it. Write the reflected `.geo` to `new/` under the same file name,
   then run `b2ag` in `new/`. Repeat the audit on the new `b2fgmtry`. Write the block option as `--block=LO:HI` with an equals
   sign: the endpoints are often negative (`-1:72`), and argparse would otherwise read them as another option.

       cp old/b2ag.dat new/b2ag.dat
       echo "'b2agbs_average_flux' '1'" >> new/b2ag.dat
       symmetrise_carre.py --block=LO1:HI1 --block=LO2:HI2 old/case.geo new/case.geo
       cd new; b2run b2ag; cd ..
       audit_b2_reflection.py new/b2fgmtry LO1:HI1 LO2:HI2

   `b2run b2ag` takes a few seconds and ends with `STOP b2ag`. In `new/` it prints `Failed opening fort.33` and still writes a
   `fort.30`: that is expected here (the mirrored mesh does not exist yet), and **that `fort.30` is not the one to use**; step 7 makes
   the final one. Only `new/b2fgmtry` matters at this stage.
   On the STEP case this took the corner mismatch from 3e-4 m to 0, the cell-volume mismatch from 4.6% to 3e-14 and the
   `hx` mismatch from 6.8% to 1e-13.

3. Build the mirrored Eirene mesh (`new/mesh/fort.33`-`fort.35` and `new/mesh/triangle_reflection.txt`) and audit it:

       mirror_eirene_mesh.py --original old --old-b2 old/b2fgmtry --new-b2 new/b2fgmtry --output new/mesh
       audit_mesh_symmetry.py --b2-geometry new/b2fgmtry --eirene-vertices new/mesh/fort.33 --eirene-triangles new/mesh/fort.34 --output new/mesh_audit.report

4. Mirror the additional surfaces and surface sources of the formatted Eirene input. First find the mirror partner of every
   additional surface of block 3b (numbered 1..N in the file):

       make_surface_pairing.py old/fort.1 pairing.txt

   Two surfaces are partners when the end points of one are the mirror images (Z to -Z) of the end points of the other, in
   reversed order, to within 0.01 cm (`--tolerance`); a surface that is its own mirror image is paired with itself. The
   script writes nothing and lists the surfaces concerned if any surface has no partner or several, because those need a
   decision by hand. On the STEP input it finds 109 pairs and 4 self-mirror surfaces, the same as a pairing made by hand.
   `pairing.txt` has one line `i j` per surface (`i i` for a self-mirror one); lines starting with `#` are comments, so a
   hand-made file works too. The mapping must be a complete involution. Then:

       symmetrise_eirene_inputs.py old/fort.1 pairing.txt new/fort.1

   For every pair, the surface in the upper half receives exactly the mirror image of the end points of its partner, and the same
   is done for the source coordinates of the strata with a surface source; all other parameters of paired surfaces and sources
   must agree, or the script stops. The changed lines are listed. On the original STEP input this changed one additional-surface
   end point and the matching gas-puff source line by 0.001 cm, the only departures from mirror symmetry in the file.

5. Project the restart state and the neutral background onto the symmetric grid (the neutral tool reads
   `new/mesh/triangle_reflection.txt`):

       symmetrise_state.py --source old/b2fstate --old-b2 old/b2fgmtry --new-b2 new/b2fgmtry --output new/b2fstati
       symmetrise_neutral_restart.py old/fort.13 new/mesh new/fort.13

   (`--source` is the name of your starting state in `old/`: `old/b2fstati` if that is what you have; `No such file or directory`
   means the name is wrong. The state projection takes about 20 s.)

   On the STEP case the species inventories were conserved to 1e-16 and the electron and ion energies exactly.

   `symmetrise_state.py` writes `b2fstati` in the native layout (6 values of 22 characters, 13 decimals; reading and rewriting a
   native file reproduces it byte for byte). This matters: `b2run` judges a run by comparing the size of the new `b2fstate`
   with that of `b2fstati`, and reports "Run apparently failed: Optional files not moved" if the new one is smaller.

   The neutral tool reads the short-format `fort.13` of the SOLPS-ITER coupling (one record of 9 virtual species x cells x 5
   moments) and writes a moments-only restart; the mesh-clipped run needs a matching Eirene that reads it (not part of this step
   of the branch). If `old/fort.13` was written on the ORIGINAL mesh, its cell count differs from the new mesh, and the tool
   then uses `new/mesh/triangle_origin.txt` (written by `mirror_eirene_mesh.py`: the original triangle each new triangle comes
   from) to carry the cells over: the lower half keeps its values and the upper half takes the mirror image of the lower half,
   so the values of the original upper half are not used. Rerun `mirror_eirene_mesh.py` if `triangle_origin.txt` is missing.

   **On the public stepspeedup branch do not use `new/fort.13` as a restart:** it needs an Eirene that reads the saved virtual
   background in the short format, which is not in the Eirene of this branch. Start without a neutral restart (no `fort.13` in the run
   directory); the first Eirene call then starts cold and the first tens of steps are a strong transient.

6. Check the result and the boundary-condition input (each tool writes a report to the file named last):

       audit_state_symmetry.py new/b2fstati new/b2fgmtry new/state_audit.report
       measure_plasma_asymmetry.py new/b2fstati new/b2fgmtry new/asymmetry.report
       audit_boundary_symmetry.py old/b2.boundary.parameters new/boundary_audit.report

7. Make the symmetric case runnable. `b2run b2mn` rebuilds any base file (`b2fgmtry`, `fort.30`, `b2fpardf`, `b2fstati`,
   `b2frates`) that is **older than the files it is made from** (see `runs/Makefile`), and if it cannot it stops with
   `No rule to make target 'b2ai.dat'` or similar, or, worse, it replaces your symmetric state by a flat one. So do the following **in
   this order**. The base files live in `runs/<your case>/baserun`; the run directory is a sibling, `runs/<your case>/<run>`.

   a. Put the symmetric files into `baserun` (you are in `baserun` already; `new/` is a subdirectory of it):

          cp new/b2ag.dat b2ag.dat
          cp new/<name of your .geo file> .
          cp new/mesh/fort.33 new/mesh/fort.34 new/mesh/fort.35 .

      Keep your original `b2ah.dat`, `b2ai.dat`, `b2ar.dat` and `b2frates` there. Do not copy `new/b2fgmtry` or `new/fort.30`; they
      are regenerated next.

   b. Regenerate the grid files from the symmetric inputs, in this order, each in `baserun`:

          b2run b2ag      # makes b2fgmtry and fort.30. It must say "Writing fort.30 ... to match ../fort.33 format"; if it says
                          # "Failed opening fort.33" the mirrored fort.33-35 are not in baserun
          b2run b2ah      # makes b2fpardf from b2ah.dat and the new b2fgmtry
          b2run b2ai      # makes a flat b2fstati in baserun (replacing any there; your own state is in old/ and new/); it is not your
                          # starting state, but make needs the file to exist

   c. `b2frates` (the atomic rates, made from ADAS data) is now older than the new `b2fpardf`, so `b2run b2mn` would try to remake
      it with `b2ar`. Either keep your original `b2frates` (right if the species set is unchanged) and let step d mark it as current,
      or remake it with `b2run b2ar`, which needs the ADAS data in `modules/adas` of the SOLPS-ITER tree.

   d. Put the timestamps of the base files in dependency order with the stock script, run from the directory that **contains**
      `baserun` (it does `cd baserun`):

          cd ..
          correct_baserun_timestamps

   e. In the run directory (`runs/<your case>/<run>`) put the usual inputs (`b2mn.dat`, the `b2*.parameters` files,
      `user_data.input`, the JSON input files if your original run directory has them), the mirrored `fort.1` (`cp
      ../baserun/new/fort.1 fort.1`) and, **last**, the symmetric state:

          cp ../baserun/new/b2fstati b2fstati

      It must be newer than every base file of step d (if you copy it earlier, the timestamp script makes it older: copy it after).

   f. Check before submitting. This lists what `b2run b2mn` would rebuild; there must be no `b2ag.exe`, `b2ah.exe`, `b2ai.exe` or
      `b2ar.exe` line:

          gmake -n -f $SOLPSTOP/runs/Makefile b2mn.prt | grep "b2a[ghir].exe"

      A `b2ai.exe` line means the run-directory `b2fstati` is older than a base file (copy it again); a `b2ar.exe` line means
      `b2frates` is older than `b2fpardf` (repeat step d).

   The surface and source physics are not certified by these tools: audit them separately. How to submit and which settings to use
   are in the user guide of the branch (`stepspeedup_user_guide.md`, section 4).

## Settings added by this branch (b2ag.dat and b2mn.dat)

| Key | Default | Effect |
|---|---|---|
| `b2agbs_average_flux` (in `b2ag.dat`) | 0 | 1 takes the reference poloidal flux of each flux tube as the average over all valid western faces of the segment, so that the corrected `pbsx`, `hx`, `hy` do not depend on which face is chosen and are invariant under the up-down reflection (B2.5 branch `feature/moultod/stepspeedup`). Only active if `b2agdr_redef_gmtry` is not 0 (the default is 1). For connected double null grids made symmetric with steps 0-3. Documented in `b2input.xml` and `doc/solps/solps.tex`. |
| `eirene_symmetrise_sources` | 0 | 1 averages the sources that Eirene hands to B2.5 with their up-down mirror image (B2.5 branch `feature/moultod/stepspeedup`). **Only meaningful for a connected double null** on an up-down symmetric grid; the key is read at the first Eirene call, before any Eirene run, and the run stops if the grid is not a connected double null (two cuts, same `topcut` and `bottomcut`, two poloidal blocks). Whether the grid itself is symmetric is not checked: use steps 0-3 above. Particle, energy and linearisation coefficients are even; the constant part of the parallel-momentum source is odd. The averaging is applied once per Eirene call, to the final sources, and prints `EIRENE_SOURCES_SYMMETRISED` in the log. With 0 the code is unchanged. Documented in `b2input.xml` and `doc/solps/solps.tex`. |

The key is read once per run, so it cannot be switched during a run. It averages the sources of each cell with those of its
mirror cell. It does not replace the symmetric grid, mesh and `fort.1` of steps 0-4, which the mirror cells need to be equal.


## Between-call source reconstruction and ramp (b2mn.dat)

These act on the steps without an Eirene call (`eirene_mod` > 1 in `b2.neutrals.parameters`) and need Eirene coupled
in-process. All default to off, in which case the code is unchanged. They are independent of the grid symmetrisation.

| Key | Default | Effect |
|---|---|---|
| `eirene_rate_reconstruct` | 0 | 1: scale the ionisation, recombination and electron-cooling sources by the ratio of Eirene's own summed rates at the present and last-call plasma state, per cell, and add a charge-exchange momentum and ion-energy response. |
| `eirene_rate_ratio_min`, `eirene_rate_ratio_max` | 1e-3, 1e3 | Limits of those ratios (min > 0). |
| `eirene_rate_fixed_ne` | 0 | 1: keep the electron density of the last call in the rates, so the scaling follows Te only. Not tested separately from the other settings. |
| `eirene_ramp_steps` | 0 | N > 0: blend the previous source into the new one over N steps after each call. |
| `eirene_ramp_state_file` | none | Normally not needed: with `eirene_ramp_steps` > 0 the run writes `../eirene_ramp_state` (run directory) at each checkpoint and at the end, and the next run loads it automatically if its simulation time matches the restart state; otherwise, or if there is no file, a warning is printed and the first call is not ramped. Set it to use another file, for example `'../eirene_ramp_state.0007'`; a named file must match or the run stops. |

Documented in `b2input.xml` and `doc/solps/solps.tex`. The scheme was developed for, and tested on, a STEP double-null case;
it has not been compared with a stock run on another case.

## Parallel speed-up and restart helpers (b2mn.dat)

| Key | Default | Effect |
|---|---|---|
| `balance_last_step` | 0 | 1 (with `balance_netcdf` 1): balance arrays updated and written only on the last step of a run (or when the remaining `b2mndr_cpu`/`b2mndr_elapsed` is below twice the longest step). The Eirene part of the balance is still gathered at every call. |
| `b2stbr_sna_history_file` | none | Restart file with the `SNA_CORR` controller memory `t1`: first line the number of species, second line the values. |

OpenMP and idle ranks (no keys, set at build time and by the launcher):

- Build with `source SETUP/openmp` before `gmake` (see `solps.tex`). The change in `b2news_.F`, `b2npmo.F` and `b2mod_ma28_for_7diag.F` keeps
  a separate pivot state for each species in the OpenMP solves of the 7-diagonal and momentum equations and removes the restarts of the MA28
  solver that discarded it.
- In an MPI run, ranks 1 and up wait for the next Eirene call with a non-blocking broadcast and a 1 ms sleep between polls (`b2mod_main.F`,
  `xertst.F`), so that they do not occupy cores that the OpenMP threads of rank 0 need. **This needs the matching Eirene change** (`eirsrt.F`
  in `src/interfaces/couple_SOLPS-ITER`, branch `feature/moultod/stepspeedup` of `dmoultonwork/EIRENE`): blocking and non-blocking
  broadcasts do not match in Open MPI, so build B2.5 and Eirene from these branches together.
- Thread count: only rank 0 should run several threads. With Open MPI the number given to `OMP_NUM_THREADS` is not necessarily passed to rank 0,
  so check in the log how many threads are in use. The idle-rank sleep interval is fixed at 1 ms and has been tested with Open MPI only.
- Tested on one STEP double null case: the plasma with 8 threads agreed with that of 1 thread and of the serial code. Not tested with other
  MPI libraries or cases.
