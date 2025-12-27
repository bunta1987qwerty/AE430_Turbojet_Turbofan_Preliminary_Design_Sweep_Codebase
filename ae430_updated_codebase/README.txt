AE430 Engine Design Codebase (Updated)

What's included:
- ae430/ : Python 2.7 package with turbojet + turbofan cycle analysis
- config_default.json : updated baseline config
- sweeps/ : example sweep definitions for turbojet and turbofan

Typical commands:
  python -m ae430.run_pipeline --config config_default.json
  python -m ae430.run_sweep --mode mp --nproc 8 --sweep sweeps/turbojet_sweep.json --outdir results --clean
  python -m ae430.plot_sweep_3d --csv results/summary_all.csv --x param.design.pi_c --y param.spools.omega_comp --z param.inlet.pi_d_factor --c res.fc --filter-feasible
