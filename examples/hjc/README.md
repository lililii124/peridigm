# HJC correspondence material

`HJC Correspondence` adds pressure-dependent strength, strain-rate sensitivity,
irreversible crushing and scalar material damage to the existing corotational
correspondence interface. The integration conventions are described below.
Density, stress and time units must be consistent.

## Homogeneous verification example

```sh
python3 generate.py --output run
cd run
/path/to/Peridigm hjc.xml
# Or: mpiexec -np 2 /path/to/Peridigm hjc.xml
```

The generator uses only the Python standard library. It creates a 20 mm cube
of 8 x 8 x 8 points, 2.5 mm spacing, with horizon 3.01 times the spacing.
The fixed timestep is 0.1 microseconds and the duration is 100 microseconds.
All points follow the same affine compression/shear cycle:

```text
h = sin(pi*t/100 microseconds)^2
F = [1-.015h    .04h       0    ]
    [   0     1-.015h    .02h  ]
    [   0        0     1-.015h ]
```

This is a homogeneous material/kinematic patch test. It exercises compaction,
plastic shear, damage and unloading, without free-surface dynamics or crack
localization. The supplied parameters demonstrate the model; they are not a
calibration. The Exodus output includes `HJC_Pressure` (compression positive),
`Von_Mises_Stress`, `HJC_Damage`, `HJC_Plastic_Volume`,
`Equivalent_Plastic_Strain`, `HJC_Density_Measure` and `HJC_Acoustic_Modulus`.

![Pressure, equivalent stress and damage histories](response.png)

The material-point reference uses an independently checked constitutive update,
independent polar/log-strain kinematics and a 5 ns timestep.
Pressure/stress peak-normalized RMSE was 0.0248%/0.0393% at 100 ns and
0.0117%/0.0241% at 50 ns. Serial and two-process histories differed by less
than 1e-6 Pa in pressure and equivalent stress. These are patch-test results,
not experimental validation of concrete parameters.

Use separate output directories for `--dt 5e-8`, `--cells 12`, another
`--horizon-ratio`, or `--hourglass`. The last option controls the existing
correspondence stabilization. A homogeneous patch does not establish the
stability of nonuniform correspondence solutions or mesh-objective softening.
The base material integrates a corotational rate; rapidly rotating prescribed
motion also requires timestep refinement to control kinematic truncation error.
Peridigm's initial critical-step estimate uses the initial bulk modulus;
the cubic EOS can become stiffer during compression. Specify a sufficiently
small `Fixed dt` and check timestep refinement for dynamic applications.

## Sphere impact on an HJC disk

`disk_impact.yaml` adapts the repository's existing
[`disk_impact`](../disk_impact/disk_impact.yaml) example. It reuses the same
`disk_impact.g` mesh: a 74 mm diameter, 2.5 mm thick disk and a 10 mm diameter
elastic sphere (14,445 disk points and 2,522 sphere points). The sphere starts
at -100 m/s along z. Both bodies are otherwise free; no motion is prescribed
after the initial conditions.

The disk uses the HJC demonstration parameters above. Its original critical
stretch bond-damage model is removed, so the example isolates HJC material
damage without an additional bond-breaking law. The elastic sphere, horizons,
contact radius and contact stiffness follow the upstream example. An initial
rigid translation moves the sphere 8 mm closer to the disk, leaving a 0.5 mm
surface gap and shortening the free-flight interval. Contact search runs every
10 steps; a fixed 50 ns timestep and a 50 microsecond end time resolve the early
impact, wave propagation and unloading response.

Run from this directory using a YAML-enabled build:

```sh
/path/to/Peridigm disk_impact.yaml
python3 plot_impact.py hjc_disk_impact.e --labels "50 ns"
```

![Impact force, sphere velocity and disk damage](impact_response.png)

The plotter needs NumPy, netCDF4 and Matplotlib. It also writes CSV histories
and a JSON summary, checks positive disk deformation determinants, finite
fields, irreversible histories and free-body momentum balance. `Contact_Force`
already includes point volume, so the sphere's plotted z-force is the direct
sum of its z-components. Momentum
is computed as a vector sum from point velocities and masses; the existing
`Global_Linear_Momentum` diagnostic sums block-wise magnitudes and is not used
for this conservation check.

For a timestep comparison, copy the YAML, set `Fixed dt` to `2.5e-8`, change
`Output Frequency` to `40` and use a different output filename. Then pass both
Exodus files to `plot_impact.py --labels "50 ns" "25 ns"`. Keeping the search
frequency at 10 steps also halves its physical search interval. For MPI, first
decompose the existing mesh with the usual SEACAS `decomp` command; use `epu`
to merge partitioned output before plotting.

On this mesh, the 50 ns/25 ns histories have peak-normalized RMSE of 0.416%
in contact force and 0.391% in volume-averaged disk damage. The final sphere
velocities differ by 0.0412 m/s. The minimum disk deformation determinant
over both runs is 0.9129, and relative vector-momentum drift is below 4e-15.
Serial/four-process runs over the first 20 microseconds agree to 1.2e-10 N
in saved contact forces. These checks establish numerical consistency for
this example; they do not establish experimental accuracy.

### Impact contours

```sh
python3 plot_contours.py hjc_disk_impact.e
```

![HJC damage and equivalent stress during sphere impact](impact_contours.png)

The figure shows the central impact region at 10, 25 and 50 microseconds from
the 25 ns run. Each field uses one fixed color scale across all frames. The
reference y=0 section retains the y>=0 half of the disk beneath the complete,
neutral-colored elastic sphere. The orthographic view uses a
displacement scale of one and is cropped to the impact region.

Material-point values are mapped by `Element_Id` to the original mesh cells
without smoothing the fields. For visualization only, displacements are
volume-averaged at shared mesh nodes and interpolated along cut edges. Fully
damaged cells remain visible. The plotter uses the same dependencies as the
response plotter; `--times` selects saved times in microseconds, `--mesh`
selects the original mesh, `--overview` shows the entire half disk, and
`--output` sets the image path.

At 50 microseconds, within 15 mm of the impact axis, timestep refinement gives
volume-weighted RMS differences of 0.00436 in damage and 1.40 MPa in equivalent
stress. Maximum pointwise differences are 0.0960 and 47.5 MPa, respectively;
the small global-history differences do not establish convergence of local
softening fields.

This is an early-time impact demonstration, not a perforation or fragmentation
benchmark. HJC material damage does not remove particles or break bonds. The
demonstration parameters are not calibrated to the original brittle disk.

## Parameters and history

The generated XML lists all required parameters. `Shear Modulus`,
`Compressive Strength`, `Tensile Strength`, `Crushing Pressure`,
`Locking Pressure`, `K1`, `K2` and `K3` have stress units.
`Reference Strain Rate` has inverse-time units. `A, B, C, N, EFMIN, SFMAX,
D1, D2`, both volumetric strains and `Hourglass Coefficient` are dimensionless.
The initial bulk modulus is `Crushing Pressure / Crushing Volumetric Strain`;
an optional `Bulk Modulus` must agree. `Locking Plastic Volumetric Strain`
is the HJC `UL` plastic offset, not the total strain at full compaction.

The six two-step history fields preserve log volume, irreversible plastic
volume, maximum compression, density measure, damage and equivalent plastic
strain. Trial evaluations always start from `STEP_N`, so repeated residual
evaluations do not compound damage. The base correspondence material supplies
the unrotated deformation rate, rotates Cauchy stress, and assembles forces.

The integration choices are:

- `mu=exp(-integral(trace(D)*dt))-1`; the previous density measure chooses the
  EOS branch and unloading modulus. The compacted EOS is cubic in
  `(mu-UL)/(1+UL)`. Its transition is not forced to be continuous.
- The pressure cutoff uses old damage; tensile normalized strength is
  `A*(1-D-p/T)` when `p<0`.
- Rates below the reference rate retain the negative logarithm, with a
  `1e-8` floor in the chosen time units. Final strength is bounded between
  zero and `SFMAX`.
- Damage increases only if `D1*(p/FC+T/FC)^D2 > EFMIN`; `EFMIN` is an
  activation gate, not a denominator floor. When this gate is initially open,
  damage starts at `min(1,1e-4/[D1*(T/FC)^D2])`, while plastic strain is zero.
  The plastic-strain counter freezes at exactly zero strength.

These conventions define this implementation. `HJC_Damage` is material damage
and is independent of peridynamic `Bond_Damage`. This material does not delete
points or break bonds; fully damaged material retains confined compressive
strength. Thermal expansion is not supported.

## Tests and references

`ctest -R utPeridigm_HJCCorrespondenceMaterial --output-on-failure` runs
analytic invariants, engineering-shear checks, damage/crushing tests,
parameter rejection, six independent reference transitions and DataManager
history/tensor-mapping checks. The compact fixtures reconstruct the starting
state and trial deviator from numerical output; they verify one-step
constitutive returns. The fixtures are self-contained numerical data.

1. Holmquist, T. J., Johnson, G. R., and Cook, W. H. (1993).
   *A computational constitutive model for concrete subjected to large strains,
   high strain rates, and high pressures*. Proceedings of the 14th
   International Symposium on Ballistics, Quebec, pp. 591-600.
