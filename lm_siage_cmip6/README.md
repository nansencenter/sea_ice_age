# Lagrangian sea-ice age for CMIP6 model data

This directory adapts the established Lagrangian sea-ice-age calculation to
CMIP6 model output. The original observational notebooks in the repository
root are reference material and are not modified.

The current production test is NorESM2-MM historical `r1i1p1f1`, 1850–1859.
Its native-grid velocity components are rotated with the existing
perturbation/Jacobian method, remapped as Cartesian tangent vectors to the
432 × 432 Arctic LAEA grid, and then supplied to the age algorithm.  Internal
mesh distances and drift displacements are in kilometres and kilometres per
day; model drift is converted from m s-1 using `1 m s-1 = 86.4 km day-1`.

## Input preparation

Put matching CMIP6 daily files in `data/input/`.  The files may cover one or
several years, but every selected year must be contained within one file per
variable.  The three variables must have identical dates:

```text
siconc_SIday_NorESM2-MM_historical_r1i1p1f1_gn_18500101-18591231.nc
siu_SIday_NorESM2-MM_historical_r1i1p1f1_gn_18500101-18591231.nc
siv_SIday_NorESM2-MM_historical_r1i1p1f1_gn_18500101-18591231.nc
```

The filename is parsed from the right, so model names containing hyphens are
safe.  Filename fields are checked against the NetCDF global metadata.  The
source series must be unique and daily-contiguous.

Choose the dataset in
`config/noresm2_mm_historical_r1i1p1f1.yaml`:

```yaml
dataset:
  model: NorESM2-MM
  experiment: historical
  variant: r1i1p1f1
  table: SIday
  grid_label: gn
```

Only files matching all five fields are processed.  This is the intended place
for a user to change model, experiment, variant, table, or native-grid label.
Use model/experiment/variant-specific work paths when creating another config.

The same config contains paths to the target LAEA grid and the established
perturbation remapper.  Change those two absolute paths on another computer or
on NIRD.  Each native model grid needs its own weights file; the first remapped
year creates it and later years reuse it.

## Environment and checks

```bash
conda env create -f environment.yml
conda activate lm-siage-cmip6
pytest -q
```

The BAMG extension must match Python 3.11 and the platform.  The current macOS
build cannot be copied to Linux; build BAMG under Python 3.11 on NIRD.
BAMG uses randomized vertex insertion by default. Although its C++ code resets
a fixed seed, the C-library `rand()` sequence is platform-specific, so macOS
and Linux meshes are not expected to be bit-for-bit identical. Each platform's
annual restart files must therefore be continued on that platform rather than
replaced with checkpoints from another build.

List the selected source years without writing data:

```bash
python scripts/00_remap_model_data.py \
  --config config/noresm2_mm_historical_r1i1p1f1.yaml \
  --dry-run
```

## Full run

Run native-grid preprocessing and every age stage with:

```bash
python scripts/run_pipeline.py \
  --config config/noresm2_mm_historical_r1i1p1f1.yaml
```

The stages are:

```text
00_remap_model_data.py       select, rotate, and annually remap model input
00_make_mesh_model.py        build the initial triangular mesh
02_move_mesh_model.py        advect/remesh it with projected drift
03_interpolate_siconc_model.py
04_find_autumn_minima_model.py
05_compute_age_model.py
11_grid_age_model.py
12_netcdf_age_model.py       collect daily fields into complete annual files
```

The first model year is an initialization year: age fields begin just after the
first autumn-minimum cohort is established.  Subsequent complete years contain
the full daily series.

The original age calculation obtains individual age fractions by differencing
cumulative autumn-minimum cohorts.  Those cumulative fields must be nested:
older ice must be a subset of every newer cohort, and the newest cohort must not
exceed observed concentration.  The model workflow enforces that physical
invariant by capping older cumulative cohorts by the next-newer cohort before
the differences are calculated.  This does not add ice; it prevents negative
age fractions when an independently calculated annual minimum is lower than a
propagated older-ice estimate.

## Observational safeguards and model data

Model fields remove several sources of observational error, but they are exact
only at the native model grid points. Vector rotation, remapping, interpolation
to moving elements, and remeshing are still numerical approximations. The
workflow therefore separates data-quality treatments from mesh-conservation
safeguards.

Treatments that are not carried into the model workflow include the
observation-only uncertainty notebooks, satellite quality screening,
SIC-conditioned three-dimensional gap filling and Gaussian smoothing of drift,
and the hand-drawn coastal/domain exclusions used by the observational mesh.
Model land/status fields and the actual model ice domain are used instead.

The following are retained because they protect the Lagrangian calculation,
not because observations are noisy:

- conversion of east/north drift to LAEA `dx/dt, dy/dt` in km day-1;
- BAMG remeshing and conservative source-to-destination element mapping;
- area-ratio correction for convergence and divergence;
- capping transported ice by the model concentration at the destination time;
- fixed outer mesh nodes, provided the buffered boundary remains outside ice;
- local velocity halving only when a full daily step would invert a triangle
  in floating-point or BAMG's integer-coordinate geometry;
- non-negative, nested cumulative age cohorts.

### Model-data switches

The observational and revised model choices are switches rather than hidden
code changes:

```yaml
preprocess:
  drift_extent_mode: none
  remap_neighbours: 4
  complete_drift_on_ice: true
  drift_completion_concentration_threshold: 0.10
  drift_completion_max_distance_km: 50.0
  drift_completion_donor_domain: non_land
  unresolved_ice_drift_policy: exclude

variables:
  concentration: siconc_for_age
  model_concentration: siconc
  drift_support_mask: drift_supported

ice:
  sic_interpolation_mode: model_linear  # or legacy_zero_fill
  missing_sic_policy: error_on_ice_domain  # or error / zero
  autumn_search_start: "08-01"
  autumn_search_stop: "09-15"

mesh:
  drift_interpolation_mode: model_linear  # or legacy_zero_fill
  missing_drift_policy: nearest_on_ice      # or error / error_on_ice / zero
  nearest_drift_max_distance_km: 200.0

age:
  clip_autumn_area_ratio: false
  autumn_area_ratio_clip: [0.97642671, 1.02437728]
```

`model_linear` still performs the unavoidable interpolation from the fixed
LAEA cells to moving mesh nodes/elements. It does not perform an observational
coastal correction: bilinear weights are renormalized over finite model-ocean
cells, so land/fill cells are not treated as zero-concentration or zero-drift
ocean.

The remapper preserves the original remapped concentration in `siconc` and
writes the explicit `drift_supported` mask. Where concentration is at least
10%, missing drift is completed only from a non-land donor within 50 km.
Unresolved cells remain visible in `siconc` but are masked from
`siconc_for_age`; no zero drift is silently assigned. Annual remapped files
record repair, exclusion, and repair-distance diagnostics in their global
attributes.

The fully strict `error` test showed that NorESM2-MM drift is not defined over
the entire ice-free buffer mesh. On the first test day, 5,136 non-boundary mesh
nodes lacked drift, but none touched an element with SIC ≥15%. Similarly, 584
elements lacked SIC, all where maximum SIC over 1850–1859 was exactly zero.
After two days of motion, four nodes touching nearly 100% SIC moved into small
target-grid drift gaps; their nearest finite remapped drift was 25–37 km away.
Because NorESM SIC and drift use different lon/lat grids, they must not be
compared at the same native array indices. A valid audit on the common LAEA
grid found missing drift in 1,462,217 of 72,180,631 ice-bearing model-ocean
cell-times (2.03%); within nominal status cells the fraction was 446,448 of
53,446,749 (0.835%). Missing drift therefore cannot safely be assumed absent
for this remapped model input.

The recommended SIC policy therefore stops inside the potential ice domain and
sets only verified zero-ice buffer gaps to zero. For drift,
`nearest_on_ice` repairs only nodes touching ice, only from finite remapped
model drift, and only within the configured 200 km bound; it stops beyond that
distance. This is a guarded response to the demonstrated native staggered or
coastal velocity mask, not a claim that model drift is intrinsically uncertain.
Missing drift on inactive buffer nodes is zero because no ice or age is
transported there. The less permissive `error_on_ice` and fully strict `error`
switches remain available for auditing another model. Always start a new model
with `error`; enable bounded repair only after quantifying the source gaps.

The remapper now repeats the static target `status_flag` on every daily record.
The inherited remapper wrote it only at time index zero and left NetCDF fill
bytes in later records. This is a metadata/mask repair only and does not alter
SIC, drift, or the target status classes.

The area ratio represents convergence or divergence of a moving element. It is
needed in cohort transport and remains enabled there. The original autumn-only
limits clipped the ratio to approximately ±2.4%, most likely to prevent noisy
or poorly mapped observational deformation from dominating the short minimum
calculation. For internally consistent model drift, that clipping suppresses
physical deformation and makes autumn propagation inconsistent with the daily
age propagation, which uses the full ratio. The recommended model test
therefore sets `clip_autumn_area_ratio: false`; the limits remain in the config
so a direct A/B comparison is possible.

The revised minimum window is 1 August–15 September. Because all candidate
fields are propagated to the window endpoint, 15 September remains the annual
age-class reference date used by the observational product. The first age
output is therefore 16 September of the initialization year. On 15 September
the new minimum is stored as the next day's cohort seed but is not included
prematurely in that day's age. This matches the original notebook's annual
rollover and keeps at most six cumulative transported cohorts and seven
exported age fractions.

Velocity-halving events remain a separate mesh safeguard and are reported
daily. In the complete NorESM2-MM 1850–1859 switch test they were frequent,
often requiring tens of local correction iterations and sometimes affecting
hundreds to thousands of nodes. The run remained valid, but this is too common
to treat as an occasional fallback in a centuries-long production archive.
Before production, compare this daily-step result with a sub-daily advection
sensitivity or revised mesh strategy; do not silently disable the inversion
protection. The safeguard checks both floating-point triangle area and BAMG's
fixed integer-coordinate representation. The second check prevents a very thin
positive triangle from changing sign through BAMG coordinate quantization
immediately before remeshing.

The completed original-setting run is preserved. The independent revised test
uses:

```bash
python scripts/run_pipeline.py \
  --config config/noresm2_mm_historical_r1i1p1f1_model_switch_test.yaml
```

Its work files and NetCDF outputs are under
`work/NorESM2-MM/historical/r1i1p1f1_model_switch_test/` and
`output/sensitivity_tests/model_data_switches/`, so it cannot overwrite the
baseline calculation.

For a limited test, use date bounds or stage bounds:

```bash
python scripts/run_pipeline.py \
  --config config/noresm2_mm_historical_r1i1p1f1.yaml \
  --to-stage remap \
  --stop 18501231

python scripts/run_pipeline.py \
  --config config/noresm2_mm_historical_r1i1p1f1.yaml \
  --from-stage mesh \
  --stop 18501231
```

## Outputs and automatic restart

Public age products are grouped by model:

```text
output/netcdf/NorESM2-MM/
  SeaIceAge_NorESM2-MM_historical_r1i1p1f1_1850.nc
  SeaIceAge_NorESM2-MM_historical_r1i1p1f1_1851.nc
  ...
```

Each file contains the daily `sea_ice_age`, `status_flag`, age-fraction
variables, coordinates, projection, time bounds, and source provenance for one
completed model year.  Files are written through `.part.nc` and renamed only
after completion; `processing_complete = 1` marks a usable checkpoint.

On a normal rerun, `run_pipeline.py` first looks for the latest complete annual
file for the selected experiment. If 1853 is the newest complete historical
file, the Lagrangian processing restarts from its last date and continues with
1854, even though input files start in 1850.

An experiment is never ignored implicitly. For a scenario run, the runner may
fall back only to the explicitly configured parent experiment. This prevents,
for example, an `ssp126` calculation from picking up an `ssp585` checkpoint for
the same model and variant.

Restart also requires the corresponding internal mesh and cohort files below
`work/<model>/<experiment>/<variant>/`.  A gridded age NetCDF cannot reconstruct
those Lagrangian states.  If they are missing, the runner stops and reports the
needed paths.  Use `--no-auto-restart` only when deliberately rebuilding from
the beginning; use `--force` only when intentionally replacing existing work.

### Historical-to-scenario handoff (2014 to 2015)

Complete the historical calculation through its final daily record, normally
31 December 2014. Keep its internal `mesh_20141231.npz`,
`age_20141231.npz`, and initial mesh until all scenarios have been seeded. The
public gridded NetCDF alone is not a restart state.

For each scenario, copy the historical config to a scenario-specific config.
Change `dataset.experiment` and every model/experiment/variant-specific work,
remapping, verification, and input path. Then add, for example:

```yaml
restart:
  parent_config: config/noresm2_mm_historical_r1i1p1f1.yaml
  parent_experiment: historical
  parent_state_date: "20141231"
  bridge_input_file: work/NorESM2-MM/ssp126/r1i1p1f1/restart/historical_20141231_bridge.nc
```

Seed `ssp126` once:

```bash
python scripts/seed_scenario_restart.py \
  --config config/noresm2_mm_ssp126_r1i1p1f1.yaml
```

The seeding command copies four items into the scenario work tree: the initial
mesh, the exactly dated historical daily mesh, the exactly dated historical
cohort/age state, and a one-record remapped historical forcing bridge. Do not
rename 31 December 2014 state as 1 January 2015. The bridge record supplies the
31 December drift needed to advect the mesh into 1 January; the first scenario
concentration then limits the transported cohorts on 1 January.

Run the scenario normally:

```bash
python scripts/run_pipeline.py \
  --config config/noresm2_mm_ssp126_r1i1p1f1.yaml
```

The first run starts from the configured historical parent. Later reruns prefer
the latest complete `ssp126` output. Repeat the config and seed command for
`ssp245`, `ssp585`, or any other branch. Each scenario must use its own work
tree and must be seeded independently from the same historical state, never
from another scenario. Use `--force` with the seeding script only when you
intentionally want to replace an existing scenario seed.

Validate an annual product with:

```bash
python scripts/validate_output.py \
  output/netcdf/NorESM2-MM/SeaIceAge_NorESM2-MM_historical_r1i1p1f1_1850.nc
```

## Remapping verification data

The preprocessing stage writes testing-only files for the first 10 remapped
time steps. The combined file contains `siu`, `siv`, and `siconc`; the separate
SIC file makes concentration checking simpler:

```text
output/verification_data/NorESM2-MM/
  RemappedInputVerification_NorESM2-MM_historical_r1i1p1f1_first10.nc
  RemappedSICVerification_NorESM2-MM_historical_r1i1p1f1_first10.nc
```

Use this compact file to inspect drift direction/magnitude and concentration on
the LAEA grid.  It is never read by the age calculation and adds no per-year
diagnostic burden to long production runs.

Optional PNG drift snapshots remain testing-only and must be requested
explicitly:

```bash
python scripts/01_plot_drift_snapshot.py \
  --config config/noresm2_mm_historical_r1i1p1f1.yaml \
  --time-index 0 5 9
```

## Server workflow: historical followed by multiple scenarios

Use the same scientific switches for the historical run and every scenario.
`seed_scenario_restart.py` checks them and refuses a mismatched handoff.

1. Clone or copy `lm_siage_cmip6` and `bamg_remesher` to the server. Leave the
   original notebook repository unchanged.

2. Create the Python 3.11 environment and build BAMG on the server platform:

   ```bash
   conda env create -f environment.yml
   conda activate lm-siage-cmip6
   python -m pip install -e .
   python -m pip install -e ../bamg_remesher
   pytest -q
   ```

3. Put or link matching native `siconc`, `siu`, and `siv` files below the
   configured source directory. Keep experiments distinguishable by their
   CMIP6 filename fields. Do not combine variants.

4. Create one historical YAML config. Set the model, `historical`, variant,
   source directory, server target-grid path, bundled generic remapper,
   model-specific weights,
   work directory, and output directory. Use the recommended model switches
   above. Every historical work path should contain
   `<model>/historical/<variant>`. Put the high-frequency `mesh`, `sic`, `age`,
   and `grid` checkpoints on node-local or high-performance scratch, not a
   cloud-synchronised directory. Put final annual NetCDFs on durable project
   storage. For example:

   ```yaml
   preprocess:
     remapper_script: scripts/remap_cmip6_sea_ice.py
     complete_drift_on_ice: true
     drift_completion_concentration_threshold: 0.10
     drift_completion_max_distance_km: 50.0
     drift_completion_donor_domain: non_land
     unresolved_ice_drift_policy: exclude

   variables:
     concentration: siconc_for_age
     model_concentration: siconc
     drift_support_mask: drift_supported

   work_dir: /cluster/scratch/$USER/lm_siage/MODEL/historical/VARIANT
   output_dir: /cluster/projects/PROJECT/lm_siage/output/netcdf

   mesh:
     initial_mesh_file: /cluster/scratch/$USER/lm_siage/MODEL/historical/VARIANT/mesh/mesh_initial.npz
     daily_mesh_dir: /cluster/scratch/$USER/lm_siage/MODEL/historical/VARIANT/mesh/daily
   ice:
     daily_sic_dir: /cluster/scratch/$USER/lm_siage/MODEL/historical/VARIANT/sic/daily
   age:
     daily_age_dir: /cluster/scratch/$USER/lm_siage/MODEL/historical/VARIANT/age/daily
     gridded_age_dir: /cluster/scratch/$USER/lm_siage/MODEL/historical/VARIANT/grid/daily
   ```

   Replace `$USER` with its literal value if the batch system does not expand
   environment variables inside YAML. Preserve at least the final daily mesh
   and age checkpoint on durable storage whenever scratch may be purged.

5. Confirm file selection without writing:

   ```bash
   python scripts/00_remap_model_data.py \
     --config config/MODEL_historical_VARIANT.yaml \
     --dry-run
   ```

6. Run one remapped year first, then inspect both first-ten-step verification
   NetCDFs and optional drift snapshots:

   ```bash
   python scripts/run_pipeline.py \
     --config config/MODEL_historical_VARIANT.yaml \
     --to-stage remap --stop 18501231

   python scripts/01_plot_drift_snapshot.py \
     --config config/MODEL_historical_VARIANT.yaml \
     --time-index 0 5 9
   ```

   Confirm SIC range/coastal behavior, drift magnitude/direction, status mask,
   and that active-ice forcing is finite before launching centuries. Start a
   new model with fully strict `error` policies; if it fails only in a verified
   zero-ice buffer, use the bounded `nearest_on_ice` and
   `error_on_ice_domain` policies documented above rather than an unconditional
   zero-fill policy.

7. Run the complete historical calculation through 31 December 2014:

   ```bash
   python scripts/run_pipeline.py \
     --config config/MODEL_historical_VARIANT.yaml
   ```

   Validate the 2014 public NetCDF and retain the internal initial mesh,
   `mesh_20141231.npz`, and `age_20141231.npz` until every scenario is seeded.

8. For each scenario, copy the historical YAML to a new config. Change
   `dataset.experiment`, all experiment-specific remapped/work/verification
   paths, and use a scenario-specific restart bridge. Do not point two
   scenarios at the same work directory. Add:

   ```yaml
   restart:
     parent_config: config/MODEL_historical_VARIANT.yaml
     parent_experiment: historical
     parent_state_date: "20141231"
     bridge_input_file: work/MODEL/ssp126/VARIANT/restart/historical_20141231_bridge.nc
   ```

9. Seed each scenario independently from historical:

   ```bash
   python scripts/seed_scenario_restart.py \
     --config config/MODEL_ssp126_VARIANT.yaml
   ```

   Repeat for `ssp245`, `ssp585`, or other branches. Never seed one scenario
   from another scenario.

10. Run each scenario normally, interactively or through Slurm:

    ```bash
    python scripts/run_pipeline.py \
      --config config/MODEL_ssp126_VARIANT.yaml
    ```

    The first run uses the historical state and the 31 December forcing bridge
    to reach 1 January 2015. Later reruns use that scenario's own latest
    complete annual file. Partial `.part.nc` files are never restart points.

11. Validate representative annual outputs and the final year. Archive public
    NetCDFs together with the YAML configs and remapping metadata. Keep at
    least the final internal checkpoint if the calculation may be extended.

## NIRD/server quick start

The first NIRD deployment established the following reliable pattern. Use a
compute node or scheduler unless the system administrators have explicitly
permitted this workload on a login node.

1. Keep software and durable annual products under the project area and put
   high-frequency checkpoints on Datapeak (or the current project scratch
   service):

   ```bash
   export REMOTE_BASE=/projects/NS9588K/SeaIceAge
   export REMOTE_WORK=/nird/datapeak/NS9588K/SeaIceAge
   export WORKFLOW_DIR="$REMOTE_BASE/lm_siage_cmip6"
   ```

   `WORKFLOW_DIR` above matches the tested direct-copy deployment. If this
   GitHub repository is cloned as a whole, set it to
   `$REMOTE_BASE/sea_ice_age/lm_siage_cmip6` instead.

2. Load the available Miniforge module and create a project-local Python 3.11
   environment. A project-local package cache avoids filling the home quota:

   ```bash
   module load Miniforge3/24.1.2-0
   export CONDA_PKGS_DIRS="$REMOTE_BASE/software/conda-pkgs"
   cd "$WORKFLOW_DIR"
   conda env create \
     --prefix "$REMOTE_BASE/software/envs/lm-siage-cmip6" \
     --file environment.yml
   conda activate "$REMOTE_BASE/software/envs/lm-siage-cmip6"
   ```

3. Build BAMG on NIRD rather than copying the macOS binary:

   ```bash
   cd "$REMOTE_BASE/bamg_remesher"
   make -C bamg
   python setup.py build_ext --inplace
   python -m pip install --no-build-isolation -e .
   pytest -q

   cd "$WORKFLOW_DIR"
   python -m pip install --no-build-isolation -e .
   pytest -q
   ```

   The verified NIRD setup produced 7 passing BAMG tests and 29 passing
   workflow tests. Re-run both test suites after changing the environment or
   rebuilding BAMG.

4. Store server-specific YAML files outside Git (for example under
   `$REMOTE_BASE/runtime/config/`). Point their high-frequency work paths at
   `$REMOTE_WORK`, their final `output_dir` at durable project storage, and set
   the minimum window consistently to:

   ```yaml
   ice:
     autumn_search_start: "08-01"
     autumn_search_stop: "09-15"
   ```

5. Confirm the selected input catalogue before launching:

   ```bash
   cd "$WORKFLOW_DIR"
   python scripts/00_remap_model_data.py \
     --config "$REMOTE_BASE/runtime/config/MODEL_historical_VARIANT.yaml" \
     --dry-run
   ```

6. For an explicitly permitted interactive run, use GNU Screen so a broken
   SSH connection does not terminate the calculation:

   ```bash
   screen -S siage-historical
   python -u scripts/run_pipeline.py \
     --config "$REMOTE_BASE/runtime/config/MODEL_historical_VARIANT.yaml"
   ```

   Detach with `Ctrl-A`, then `D`; list sessions with `screen -ls`; reconnect
   with `screen -d -r siage-historical`. Record the NIRD login host shown in
   the prompt because a Screen session belongs to that host.

7. Long historical runs can be submitted in decade endpoints. The pipeline
   automatically resumes from the latest complete annual NetCDF; `.part.nc`
   files are deliberately ignored. If a chunk fails, correct the problem and
   rerun the same command without `--force`.

8. After the complete 2014 product exists, follow the historical-to-scenario
   procedure above: create an independent config and work directory for each
   scenario, run `seed_scenario_restart.py` against the historical
   `20141231` state, and then start the scenario pipeline. Every scenario is
   seeded directly from historical, never from another scenario.

For a scheduled run, change the account and paths in the supplied template:

```bash
sbatch --export=CONFIG=config/noresm2_mm_historical_r1i1p1f1.yaml \
  slurm/run_pipeline.slurm
```
