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
