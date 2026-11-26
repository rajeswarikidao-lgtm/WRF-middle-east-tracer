# Operational Updates on `oper-isotope-v4.6.1`

**Base commit**: [`ee63caf2`](https://github.com/rajeswarikidao-lgtm/WRF-middle-east-tracer/commit/ee63caf2) (`isotope-decay-tracer-v4.6.1`, tag `rn-stageH-pass`)  
**Operational commit**: `294577a2` (`oper-isotope-v4.6.1`)  
**Scope**: Host dust-routing alignment (`chem_opt=100`), exact double-precision release/decay timing for CFL-adaptive time stepping and warm restarts, and differentiated land vs. dynamic sea-water `iodg` dry deposition (`rn_iodg_ddep_opt=3`).

---

## 1. Host Dust Emission Routing (`chem/emissions_driver.F`)

- Cherry-picked upstream fix `02335bbe` (`USE module_aerosols_soa_vbs_het, only: soa_vbs_addemiss`) so host `MADE-SORGAM` dust emissions (`dust_opt = 13`) under `chem_opt = 100` match operational dust concentrations that drive `rn_coag_opt = 1` and `rn_upt_opt = 2`.

---

## 2. Exact `REAL*8` Time Keeping for Adaptive `dt` & Warm Restarts

- **`chem/chem_driver.F` & `chem/module_input_tracer.F`**:
  - Changed `curr_secs` argument in `set_isotope_tracer` and `rn_release` from default single-precision `REAL` (`real(curr_secs)`) to double-precision `REAL*8` (`curr_secs`).
  - Passed `curr_secs + DBLE(grid%dt)` instead of `DBLE(ktau) * DBLE(grid%dt)` into `rn_decay` (`chem/chem_driver.F`) and `rn_decay_factor` inside `decay_and_release` (`chem/module_input_tracer.F`).
- **Why**: Under `use_adaptive_time_step = .true.` (e.g., `30..75 s`) and multi-segment daily restarts (`restart = .true.`, `override_restart_timers = .true.`), `ktau * dt` drifts when `dt` changes across steps, and single-precision `real(curr_secs)` loses sub-second precision past `~4.5 h` (`t_hours = REAL(curr_secs / 3600.0D0)`). Both release windows and radioactive decay are now exact across adaptive steps and segmented restarts.

---

## 3. Land vs. Dynamic Sea-Water `iodg` Dry Deposition (`rn_iodg_ddep_opt = 3`)

- **`Registry/registry.chem`**:
  - Extended `rn_iodg_ddep_opt` documentation/modes with `3 = land/sea rc`.
  - Added two `&chem` namelist controls (`max_domains`, `rh`):
    - `rn_iodg_rc_sea` (default `-1.`, units `s m-1`): when `> 0.`, prescribes a uniform sea-water surface resistance $r_{c,\text{sea}}$; when `<= 0.`, computes $r_{c,\text{sea}}$ dynamically via the two-film air–sea exchange formulation below.
    - `rn_iodg_sea_alpha` (default `5.0`, dimensionless): chemical enhancement factor $\alpha$ for $\text{I}_2$ reaction with dissolved iodide/organics in surface seawater.
- **`chem/module_dep_simple.F` (`wesely_driver`)**:
  - When `rn_opt == 1`, `p_iodg >= param_first_scalar`, and `rn_iodg_ddep_opt == 3`:
    - Identifies water/lake grid cells robustly across both USGS and MODIS land-use classifications:
      `(iland == iswater_temp) .OR. (ivgtyp(i,j) == config_flags%iswater) .OR. (config_flags%islake > 0 .AND. ivgtyp(i,j) == config_flags%islake)` (handling MODIS `ISWATER = 17`, `ISLAKE = 21` as well as `iswater_temp = 16`).
    - Over land: sets `srfres(p_iodg) = config_flags%rn_iodg_rc` (operational default `330 s m-1`).
    - Over sea/water with `rn_iodg_rc_sea > 0.`: sets `srfres(p_iodg) = config_flags%rn_iodg_rc_sea`.
    - Over sea/water with `rn_iodg_rc_sea <= 0.`: computes two-film Wanninkhof (2014) / Sander (2023) liquid-phase exchange resistance:
      - $U_{10} = \operatorname{clamp}\left(\frac{u_*}{0.40}\ln\frac{10}{\max(z_0, 10^{-6})},\, 0.5,\, 40.0\right)\text{ m s}^{-1}$
      - $H_{\text{dim}}(T_K) = 2.84 \exp\left(4300 \left(\frac{1}{T_K} - \frac{1}{298.15}\right)\right) \cdot 0.082057 \, T_K$
      - $Sc_w(T_C) = \max\left(200.0,\, 1.37\left(2116.8 - 136.25 T_C + 4.7353 T_C^2 - 0.092307 T_C^3 + 0.0007555 T_C^4\right)\right)$
      - $k_w = \frac{(1 + 0.251 U_{10}^2)\sqrt{660 / Sc_w}}{360000}\text{ m s}^{-1}$
      - $r_{c,\text{sea}} = \operatorname{clamp}\left(\frac{1}{\alpha \, H_{\text{dim}} \, k_w},\, 10.0,\, 10^5\right)\text{ s m}^{-1}$.
- **`chem/dry_dep_driver.F`**:
  - Added `CASE (2, 3)` in the `ddvel(:,:,p_iodg)` guard so Wesely velocities computed by `wesely_driver` are preserved for `opt = 3`.
  - Added `config_flags%rn_iodg_ddep_opt == 3` to the `rn_ddep_iodg` diagnostic accumulation guard.

---

## 4. Reference Operational `&dynamics` / `&chem` Namelist

```fortran
&dynamics
 tracer_opt                          = 11,
 tracer_adv_opt                      = 2,
/

&chem
 chem_opt                            = 100,
 ne_area                             = 123,
 conv_tr_wetscav                     = 0,
 iso_conv_wetscav                    = 0,
 rn_opt                              = 1,
 rn_init_aj                          = 0.,
 rn_iod_frac_part                    = 0.2,
 rn_coag_opt                         = 1,
 rn_coag_factor                      = 1.,
 rn_upt_opt                          = 2,
 rn_upt_gamma                        = 0.001,
 rn_settle_opt                       = 0,
 rn_aer_ddep_opt                     = 1,
 rn_iodg_ddep_opt                    = 3,
 rn_iodg_vd                          = 0.,
 rn_iodg_rc                          = 330.,
 rn_iodg_rc_sea                      = -1.,
 rn_iodg_sea_alpha                   = 5.,
 rn_iodg_dvj                         = 0.096,
 rn_wdep_opt                         = 1,
 rn_wdep_scale                       = 1.,
 rn_wdep_host_opt                    = 0,
 rn_wd_test_rain                     = 0.,
/
```

---

## 5. Verification Summary (`2026-11-22` `00..02 UTC`, `10 km` Gulf `d01`, `1h + 1h` Warm Restart)

- **Donor CAMS IC Transfer**: `116` host chemistry/aerosol species copied, radionuclide fields initialized cleanly at `0.0`.
- **Restart & CFL Stability**: Segment 1 (`00..01Z`) and Segment 2 (`01..02Z`, `restart = .true.`, `override_restart_timers = .true.`, adaptive `dt = 30..75 s`) completed with `SUCCESS COMPLETE WRF` and `0` CFL warnings.
- **Physical Budgets at `02:00 UTC`**:
  - All 5 state tracers active (`iodg`, `iodaj`, `iodac`, `csaj`, `csac` finite and $> 0$).
  - Adsorption (`rn_upt_opt = 2`) active: `RN_UPT_IODJ` + `RN_UPT_IODC` sum `1.461e+05 Bq/m²`.
  - Land/sea `iodg` dry deposition (`rn_iodg_ddep_opt = 3`) active: land sum `1.696e+05`, Gulf sea sum `1.124e+06`.
  - 2-hour domain activity budget closure ratios: **`I-131 = 1.0013`**, **`Cs-137 = 1.0011`**.
