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
- The model accepts optional independent side-insulation k and thickness; the design engine does not use them yet and applies the back insulation to the edges.
