# AE430_Turbojet_Turbofan_Preliminary_Design_Sweep_Codebase
This repository contains a Python-based parameter sweep tool used for an AE 430 final project.  
It evaluates candidate turbojet and turbofan designs by varying cycle/component parameters, screening constraints, and ranking feasible designs by the project cost function.

## What it does
- Loads a base design configuration (JSON)
- Expands a sweep grid specification (JSON) into many cases
- Evaluates each case (inlet → compressor → burner → turbine/cooling → nozzle → installed thrust)
- Rejects infeasible cases using component constraints
- Computes and ranks the cost function for feasible designs
- Supports multiprocessing (ranked worker outputs + merge)

## Running a sweep (example)
```bash
python -m ae430.run_sweep --mode mp --nproc 8 \
  --sweep sweeps/turbojet_coarse.json \
  --outdir results_turbojet_coarse --clean --F 350e3

  Author: Leopoldo Gutierrez
  email: lgutierrez0984@sdsu.edu
