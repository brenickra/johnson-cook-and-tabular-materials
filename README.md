# Johnson–Cook and Tabular Material Model Generator for Abaqus

This Python tool provides a clean interface to calculate mechanical material properties based on the Johnson–Cook constitutive model or a custom tabular plasticity curve. It exports directly to the Abaqus `.inp` format, enabling easy integration into FEA workflows.

## Features

- Calculate Johnson–Cook parameters (A, B, n) from basic mechanical inputs
- Generate tabular plasticity data using a Ramberg–Osgood-based approximation
- Export both models to Abaqus-compatible `.inp` files
- Visualize true and engineering stress–strain curves
- Fully interactive command-line interface supporting multiple materials

## Requirements

- Python 3.7+
- `numpy`
- `matplotlib`
- `scipy`

Install the dependencies with:

```bash
pip install numpy matplotlib scipy
```

## Usage

Run the script:

```bash
python material_model_generator.py
```

Follow the on-screen prompts to:

1. Input material properties
2. View the generated stress–strain plots
3. Choose between Johnson–Cook or tabular plasticity models
4. Export the model directly to an Abaqus `.inp` file

## Example Output

### Johnson–Cook:

```
*MATERIAL, NAME=Inconel_718
*DENSITY
  7.85E-9,
*ELASTIC
  200000., 0.3
*PLASTIC, HARDENING=JOHNSON COOK
  606.154801, 928.212409, 0.233168, 0., 0., 0.
```

### Tabular:

```
*PLASTIC
  267.134965, 0.000000
  301.039024, 0.001989
  319.948709, 0.003098
  ...
```
