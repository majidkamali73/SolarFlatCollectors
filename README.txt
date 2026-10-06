Solar Collector Design Engine v1.3
Python 3.7.2 compatible; standard library only.

NEW IN v1.3
- Multiple selectable tube, absorber, cover, insulation and adhesive records.
- Tube geometry validation: OD > ID and wall thickness is reported.
- Back insulation material and thickness are explicit inputs.
- One or two covers can be selected.
- Cover-to-cover gap is an explicit input and is reported.
- Three legacy connection selections are available:
  1) tubes bonded below absorber plate
  2) tubes bonded above absorber plate
  3) tubes in-line with absorber plate
- Adhesive k remains a direct input, matching the BASIC program. The modern 2-D model also has an adhesive thickness parameter.

IMPORTANT MODEL STATUS
- The original BASIC source explicitly has the three connection selections and asks for adhesive thermal conductivity ka.
- The BASIC source has an insulation menu: glass wool Ki=0.05, mineral wool Ki=0.09, or user-entered Ki.
- The legacy source loops over one or two covers (mm=1..2) and uses tau^mm in the absorbed solar term. v1.3 preserves this for identical/effective covers.
- For two different cover materials, v1.3 uses an explicit engineering approximation: total transmittance is the product of the two transmittances and refractive index is averaged. This is NOT a rigorous multi-layer optical/radiative solver.
- Cover gap is stored and reported but is not yet active in the legacy top-loss equation. It should not be interpreted as a validated gap optimization variable yet.
- The 2-D thermal model applies a separate relation for each of the three connection types, following the legacy F' equations: bond resistance = adhesive thickness / (ka x tube OD) for tubes below the plate; the same bond on the fin path only, with the strip under the tube at tube-wall temperature, for tubes above the plate; no bond for in-line tubes. With uniform flow the model reproduces the legacy F' relations within 0.2% on a fine grid and about 1.3% (conservative) on the coarse design grid.
- Catalog numerical values marked example_engineering are placeholders and must be verified/replaced before engineering decisions.
- Current fluid model is water only. Density, viscosity, specific heat and conductivity are interpolated from a 0-100 C table at the mean fluid temperature.

RUN
1. Put all files in one folder.
2. In IDLE, open collector_design_gui_v1_4.py.
3. Run Module > Run Module.
4. Enter project data and press OPTIMIZE.

The command-line engine can also be imported by other programs.


V1.3 IMPORTANT INPUT CLARIFICATION
- Collector tilt angle beta is an explicit user input, in degrees from horizontal.
- The legacy BASIC correlation is retained: C = 365.9*(1 - 0.00883*beta + 0.0001298*beta^2).
- The solar irradiance input is used directly as the irradiance incident on the collector plane.
- No new solar-position or annual-radiation model is added in v1.3.
- This interpretation follows the legacy BASIC source, which directly asks for Solar radiation (W/m2) and separately asks for Angle of tilt (Degree).


v1.5: hydraulic flow design can be fixed by tube velocity, optimized over a velocity range, or entered manually.
