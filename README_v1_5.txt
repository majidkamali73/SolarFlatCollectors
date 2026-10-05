Solar Collector Design Engine v1.5

Main changes:
1. Back and side insulation are separate material selections and thickness inputs.
2. The GUI shows thermal conductivity, default thickness, and data status for each selected insulation.
3. The 2-D engineering model uses separate back resistance Ub=k_back/t_back and side/edge resistance with k_side and t_side.
4. The side-insulation formulation is an engineering extension of the later 2-D model, not an original BASIC equation.
5. The original BASIC source exposed Ki and back-insulation thickness in its later reconstructed model but did not expose independent side insulation.
6. Existing CSV catalogs remain plain comma-separated text files; Excel display is not required for the program to read them.

Python 3.7.2 / Tkinter / standard library only.

Corrections after v1.5 review:
- Cover infrared emissivity now comes from the emissivity column of data/covers.csv (default 0.88 when the model is called directly). The former 1 - R - tau estimate gave about 0.08 for glass and under-predicted top loss.
- Two covers: the model applies tau**M, so the engine now passes the per-cover (geometric mean) transmittance; the second cover is no longer applied twice.
- Independent side insulation is now active in coupled_case_fast, which is the path used by the design engine and the GUI.
- GUI: the adhesive selector has its own row and no longer covers the side-insulation selector.
