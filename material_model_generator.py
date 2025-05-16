import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

def func1(m, sigy, epsy, sigu):
    """
    Auxiliary function used to solve for m in the Johnson-Cook parameter estimation.

    Parameters
    ----------
    m : float
        Intermediate parameter.
    sigy : float
        Yield strength [MPa].
    epsy : float
        Yield strain.
    sigu : float
        Ultimate tensile strength [MPa].

    Returns
    -------
    float
        Result of the nonlinear equation for solving m.
    """
    return sigy * (1.0 + epsy) - sigu * np.exp(m) * ((np.log(1.0 + epsy) / m) ** m)

def func2(n, m, sigy, epsy, sigu, ymod):
    """
    Auxiliary function used to solve for n in the Johnson-Cook parameter estimation.

    Parameters
    ----------
    n : float
        Strain hardening exponent.
    m : float
        Intermediate parameter.
    sigy : float
        Yield strength [MPa].
    epsy : float
        Yield strain.
    sigu : float
        Ultimate tensile strength [MPa].
    ymod : float
        Young's modulus [MPa].

    Returns
    -------
    float
        Result of the nonlinear equation for solving n.
    """
    term1 = sigu * np.exp(m)
    term2 = sigy * (1.0 + epsy)
    term3 = (np.power(m - sigu * np.exp(m) / ymod, n) -
             np.power(np.log(1.0 + epsy) - sigy * (1.0 + epsy) / ymod, n))
    denom = (n * np.power(m - sigu * np.exp(m) / ymod, n - 1.0) *
             (np.exp(-m) - sigu / ymod))
    return term1 - term2 - sigu * term3 / denom

def calcular_coeficientes(ymod, sigy, sigu, offset):
    """
    Calculate Johnson-Cook material model parameters A, B, n.

    Parameters
    ----------
    ymod : float
        Young's modulus [MPa].
    sigy : float
        Yield strength [MPa].
    sigu : float
        Ultimate tensile strength [MPa].
    offset : float
        Offset strain (e.g., 0.002 for 0.2%).

    Returns
    -------
    A : float
        Johnson-Cook parameter A [MPa].
    B : float
        Johnson-Cook parameter B [MPa].
    n : float
        Johnson-Cook strain hardening exponent.
    epsu : float
        Ultimate true strain.
    """
    epsy = offset + sigy / ymod

    sol_m = root_scalar(lambda m: func1(m, sigy, epsy, sigu), bracket=[0.01, 1.0], method='brentq')
    m = sol_m.root
    epsu = np.exp(m) - 1.0

    sol_n = root_scalar(lambda n: func2(n, m, sigy, epsy, sigu, ymod), bracket=[0.01, 1.0], method='brentq')
    n = sol_n.root

    term_b = (np.log(1.0 + epsu) - sigu * (1.0 + epsu) / ymod)
    if term_b <= 0:
        raise ValueError("Invalid value for logarithm in B calculation.")
    B = sigu / (n * (term_b ** (n - 1.0)) * (1.0 / (1.0 + epsu) - sigu / ymod))
    A = sigy * (1.0 + epsy) - B * (np.log(1.0 + epsy) - sigy * (1.0 + epsy) / ymod) ** n

    return A, B, n, epsu

def plotar(material, A, B, n, sigy, sigu, ymod, offset, epsu):
    """
    Plot the true and engineering stress–strain curves using the Johnson-Cook model.

    This function generates a plot comparing the engineering and true stress–strain
    responses of a material described by the Johnson-Cook model. It uses the model
    parameters A, B, n to calculate the true stress from transformed strain values,
    and derives the engineering stress accordingly. The yield and ultimate strengths
    are displayed as horizontal reference lines with annotations, and the final point
    of the true curve is marked with a labeled crosshair.

    Parameters
    ----------
    material : str
        Name of the material to be shown in the plot title.
    A : float
        Johnson-Cook parameter A [MPa].
    B : float
        Johnson-Cook parameter B [MPa].
    n : float
        Johnson-Cook strain hardening exponent.
    sigy : float
        Yield strength [MPa].
    sigu : float
        Ultimate tensile strength [MPa].
    ymod : float
        Young's modulus [MPa].
    offset : float
        Yield offset strain (e.g., 0.002 for 0.2% offset).
    epsu : float
        True strain at rupture (maximum true strain).

    Returns
    -------
    strains_eng : ndarray
        Array of engineering strain values.
    stresses_eng : ndarray
        Array of engineering stress values computed from the Johnson-Cook true stress.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    strains_eng_full = np.linspace(0, epsu, 1000)
    log_valid = np.log(1 + strains_eng_full) - sigy * (1 + strains_eng_full) / ymod

    mask_valid = log_valid > 0
    strains_eng_valid = strains_eng_full[mask_valid]
    log_valid_valid = log_valid[mask_valid]

    stresses_true_valid = A + B * np.power(log_valid_valid, n)
    stresses_eng_valid = stresses_true_valid / (1 + strains_eng_valid)
    strains_true_valid = np.log(1 + strains_eng_valid)

    strains_eng = np.insert(strains_eng_valid, 0, 0.0)
    stresses_eng = np.insert(stresses_eng_valid, 0, 0.0)
    strains_true = np.insert(strains_true_valid, 0, 0.0)
    stresses_true = np.insert(stresses_true_valid, 0, 0.0)

    strain_u = strains_true[-1]
    stress_u = stresses_true[-1]

    plt.figure(figsize=(12, 6))
    plt.plot(strains_eng, stresses_eng, 'k--', linewidth=1.2, label='Stress–Strain Engineering')
    plt.plot(strains_true, stresses_true, 'k-', linewidth=1.5, label='Stress–Strain True')

    plt.axhline(y=sigy, color='gray', linestyle='--', linewidth=1, label='Yield Strength')
    plt.axhline(y=sigu, color='gray', linestyle=':', linewidth=1, label='Ultimate Strength')
    plt.text(epsu * 0.99, sigy + 5, f'{sigy:.1f} MPa', ha='right', va='bottom', fontsize=9, color='gray')
    plt.text(epsu * 0.99, sigu + 5, f'{sigu:.1f} MPa', ha='right', va='bottom', fontsize=9, color='gray')

    plt.axvline(x=strain_u, color='black', linestyle=':', linewidth=1)
    plt.axhline(y=stress_u, color='black', linestyle=':', linewidth=1)
    plt.text(strain_u + 0.002, stress_u, f"εu = {strain_u:.4f}\nσu = {stress_u:.1f} MPa",
             va='center', ha='left', fontsize=9,
             bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", lw=0.5))

    plt.title(f'Material: {material}', fontsize=12)
    plt.xlabel('Strain')
    plt.ylabel('Stress [MPa]')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return strains_eng, stresses_eng

def exportar_inp_jc(nome, material, density, ymod, nu, A, B, n):
    """
    Export material model to Abaqus .inp format using Johnson-Cook formulation.

    Parameters
    ----------
    nome : str
        Output file name.
    material : str
        Material name.
    density : float
        Density [tonn/mm³].
    ymod : float
        Young's modulus [MPa].
    nu : float
        Poisson's ratio.
    A : float
        Johnson-Cook A parameter [MPa].
    B : float
        Johnson-Cook B parameter [MPa].
    n : float
        Johnson-Cook n parameter.
    """
    with open(f"{nome}.inp", "w") as f:
        f.write(f"*MATERIAL, NAME={material}\n")
        f.write("*DENSITY\n")
        f.write(f"  {density},\n")
        f.write("*ELASTIC\n")
        f.write(f"  {ymod}, {nu}\n")
        f.write("*PLASTIC, HARDENING=JOHNSON COOK\n")
        f.write(f"  {A:.10f}, {B:.10f}, {n:.10f}, 0., 0., 0.\n")

def exportar_inp_tabular(nome, material, density, ymod, nu, sigy, sigu, offset, epsu, npts=15):
    """
    Export material model to Abaqus .inp format using tabular stress–strain data.

    Parameters
    ----------
    nome : str
        Output file name.
    material : str
        Material name.
    density : float
        Density [tonn/mm³].
    ymod : float
        Young's modulus [MPa].
    nu : float
        Poisson's ratio.
    sigy : float
        Yield strength [MPa].
    sigu : float
        Ultimate strength [MPa].
    offset : float
        Offset strain.
    epsu : float
        Strain at rupture.
    npts : int, optional
        Number of tabular points, by default 15.
    """
    eps_yield = offset + sigy / ymod
    eps_rupt = epsu
    eps_yield_pl = eps_yield - sigy / ymod
    eps_rupt_pl = eps_rupt - sigu / ymod
    n = (np.log10(sigu) - np.log10(sigy)) / (np.log10(eps_rupt_pl) - np.log10(eps_yield_pl))
    logK = np.log10(sigu) - n * np.log10(eps_rupt_pl)
    K = 10 ** logK

    sigma_nom = np.linspace(sigy, sigu, npts - 1)
    strain_nom = sigma_nom / ymod + (sigma_nom / K) ** (1 / n)
    strain_true = np.log(1 + strain_nom)
    stress_true = sigma_nom * (1 + strain_nom)
    strain_plastic = strain_true - (stress_true / ymod)

    # First special point: extrapolated stress at strain = sigy / E
    strain_true_first = sigy / ymod
    strain_plastic_first = 0.0
    s1, s2 = stress_true[0], stress_true[1]
    e1, e2 = strain_true[0], strain_true[1]
    stress_true_first = s1 - ((s2 - s1) / (e2 - e1)) * (e1 - strain_true_first)

    stress_true = np.insert(stress_true, 0, stress_true_first)
    strain_plastic = np.insert(strain_plastic, 0, strain_plastic_first)

    with open(f"{nome}.inp", "w") as f:
        f.write(f"*MATERIAL, NAME={material}\n")
        f.write("*DENSITY\n")
        f.write(f"  {density},\n")
        f.write("*ELASTIC\n")
        f.write(f"  {ymod}, {nu}\n")
        f.write("*PLASTIC\n")
        for s, e in zip(stress_true, strain_plastic):
            f.write(f"  {s:.10f}, {e:.10f}\n")

def main():
    while True:
        print("\n--- Material Model Input ---")
        material = input("Material name: ").strip()
        ymod = float(input("Young's modulus [MPa]: "))
        sigy = float(input("Yield strength [MPa]: "))
        sigu = float(input("Ultimate strength [MPa]: "))
        offset = float(input("Offset strain (e.g., 0.002 for 0.2%): "))
        epsu = float(input("Strain at rupture (e.g., 0.2 for 20%): "))
        nu = float(input("Poisson's ratio: "))
        density = float(input("Density [tonn/mm³] (e.g., 7.85e-9): "))

        try:
            A, B, n, _ = calcular_coeficientes(ymod, sigy, sigu, offset)
        except Exception as e:
            print(f"Error in calculations: {e}")
            continue

        print(f"\nJohnson-Cook parameters for {material}:")
        print(f"A = {A:.6f} MPa")
        print(f"B = {B:.6f} MPa")
        print(f"n = {n:.10f}")

        strains, stresses = plotar(material, A, B, n, sigy, sigu, ymod, offset, epsu)

        export = input("Do you want to export the material to a .inp file? (y/n): ").strip().lower()
        if export == 'y':
            opcao = input("Choose model type (1 = Johnson-Cook, 2 = Tabular): ").strip()
            nome_arquivo = material.replace(" ", "_") + ("_JC" if opcao == '1' else "_TABULAR")
            if opcao == '1':
                exportar_inp_jc(nome_arquivo, material, density, ymod, nu, A, B, n)
                print(f"File '{nome_arquivo}.inp' exported successfully.")
            elif opcao == '2':
                exportar_inp_tabular(
                    nome_arquivo,
                    material,
                    density,
                    ymod,
                    nu,
                    sigy,
                    sigu,
                    offset,
                    epsu
                )
                print(f"File '{nome_arquivo}.inp' exported successfully.")
            else:
                print("Invalid option. No file exported.")

        cont = input("Do you want to enter another material? (y/n): ").strip().lower()
        if cont != 'y':
            print("Exiting.")
            break


if __name__ == "__main__":
    main()
