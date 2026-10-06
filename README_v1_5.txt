Solar Collector Design Engine v1.5

New in v1.5:
- Flow is no longer forced to be a raw user input.
- fixed_speed: derives total mass flow from a selected mean tube velocity and the actual tube ID and number of parallel tubes.
- optimize_speed: scans a user-defined velocity range and includes velocity as a design variable, without requiring a separate mass-flow guess.
- manual_flow: preserves expert direct mass-flow entry.
- Results report mass flow, L/min, mean tube velocity, and mean Reynolds number.
- The velocity is a hydraulic design variable/check, not a claim of one universal optimum.
- Default scan 0.20 to 0.80 m/s in 0.10 m/s steps is a configurable starting range. Verify project-specific limits.

The legacy BASIC correlations remain in the thermal model. Engineering hydraulic/network extensions are reported separately.

Corrections after review:
- Cover infrared emissivity now comes from the emissivity column of data/covers.csv (default 0.88 when the model is called directly). The former 1 - R - tau estimate gave about 0.08 for glass and under-predicted top loss.
- Two covers: the model applies tau**M, so the engine now passes the per-cover (geometric mean) transmittance; the second cover is no longer applied twice.
- Side insulation material and thickness are separate inputs in the engine and the GUI. This is an engineering extension, not an original BASIC input; when not given, the back insulation is used for the edges.
- Edge loss keeps the legacy relation with (L1+L2), as in the original BASIC program: full perimeter 2*(L1+L2) with the mean temperature drop across the side insulation taken as half of plate minus ambient. The edge height L3 is set by the back-insulation thickness, and the side-insulation thickness appears only in the conduction path.
- Speed: the plate equations are solved line by line (exact tridiagonal solve across the tubes), header losses use running sums, the collector length is found by false position, and layouts wider than the limit are skipped before solving. Results are unchanged within the solver tolerance.
- GUI: the calculation runs in a background thread, so the window stays responsive; the status bar shows progress and a Stop button cancels the run.
- Tube-to-plate connection: the three legacy arrangements now each have their own relation in the model, and the bond width is the tube outside diameter (formerly the full circumference, one relation for types 1 and 2). With ordinary adhesive (k = 0.2 W/m.K) this lowers the predicted useful heat by about 11%.
