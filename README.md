# Johnson–Cook and Tabular Material Model Generator for Abaqus

This Python tool generates material input cards for Abaqus based on the Johnson–Cook constitutive model or a power-law tabular plasticity approximation. It supports both interactive and batch processing modes, exports `.inp` files automatically, and provides clean, standardized plots of stress–strain behavior.

## Features

* Compute Johnson–Cook parameters (A, B, n) from basic inputs
* Estimate tabular plasticity curves using Ramberg–Osgood-style fitting
* Export both models in Abaqus-compatible `.inp` format
* Generate consistent plots (true & engineering stress–strain)
* Choose between:

  * Interactive CLI mode for one material at a time
  * Batch mode using a CSV file (`material_data.csv`) for multiple materials

## Requirements

* Python 3.7+
* `numpy`, `matplotlib`, `scipy`, `pandas`

Install the dependencies:

```bash
pip install numpy matplotlib scipy pandas
```

## Usage

### Interactive mode

Run the script:

```bash
python material_model_generator.py
```

Follow the prompts to input material data, view plots, and export Abaqus `.inp` files.

### Batch mode

Create a file called `material_data.csv` in the script folder with the following format:

```csv
Materialname;E;Sy;Su;offset;epsolon_u;poisson;density
Steel AISI 1045;210000;450;600;0.002;0.15;0.3;7.85e-9
Aluminum 6061-T6;70000;275;310;0.002;0.12;0.33;2.70e-9
Titanium Ti-6Al-4V;113000;830;900;0.002;0.10;0.34;4.43e-9
```

Then run the script and select batch mode. It will automatically:

* Calculate Johnson–Cook and tabular models
* Export `.inp` files
* Save the plots as PNG files (no visualizer opened)

## Output Example

### Johnson–Cook (.inp)

```
*MATERIAL, NAME=Steel_AISI_1045
*DENSITY
  7.85e-9,
*ELASTIC
  210000., 0.3
*PLASTIC, HARDENING=JOHNSON COOK
  500.0, 600.0, 0.22, 0., 0., 0.
```

### Tabular (.inp)

```
*PLASTIC
  270.0, 0.000000
  289.1, 0.001112
  307.8, 0.002334
  ...
```

## Notes

* All images and material cards are saved in the same directory as the script
* Plot files are named using the pattern `<material_name>_JC.png` and `<material_name>_TABULAR.png`

## License

MIT License. Free to use and modify.
