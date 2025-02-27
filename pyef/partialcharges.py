import numpy as np
from pyscf import gto, tools
from pyscf.tools import molden

class IterativeHirshfeld:
    def __init__(self, molden_file, rad_file):
        """
        Initialize the IterativeHirshfeld class.
        Parameters:
            molden_file: Path to the .molden file containing molecular data.
            rad_file: Path to the .rad file containing pro-molecular densities.
        """
        self.molden_file = molden_file
        self.rad_file = rad_file
        self.mol = None
        self.mo_coeff = None
        self.occ = None
        self.promolecular_densities = None

    def load_molden_data(self):
        """
        Load molecular data from the .molden file.
        """
        self.mol = gto.Mole()
        self.mol.build(parse_arg=False)
        # Read molecular geometry from the .molden file
        with open(self.molden_file, "r") as f:
            molden.read(self.mol, f)
            # Rebuild the molecule after reading from the .molden file
            self.mol.build()

        # Load molecular orbital coefficients and occupations
        with open(self.molden_file, 'r') as f:
            molden_data = molden.read(self.mol, f)
            self.mo_coeff, self.occ = molden_data

    def load_promolecular_densities(self):
        """
        Load pro-molecular densities from the .rad file.
        """
        atomic_densities = {}
        with open(self.rad_file, "r") as f:
            current_atom = None
            current_density = []
            for line in f:
                if line.startswith("# Atom:"):
                    if current_atom and current_density:
                        atomic_densities[current_atom] = np.array(current_density)
                    current_atom = line.split(":")[1].strip()
                    current_density = []
                elif line.strip():
                    current_density.append(float(line.strip()))
            if current_atom and current_density:
                atomic_densities[current_atom] = np.array(current_density)

        # Match atomic densities to molecule's atoms
        self.promolecular_densities = []
        for i in range(self.mol.natm):
            symbol = self.mol.atom_symbol(i)
            if symbol not in atomic_densities:
                raise ValueError(f"No density found for atom {symbol} in {self.rad_file}")
            self.promolecular_densities.append(atomic_densities[symbol])

    def calculate_density_matrix(self):
        """
        Calculate the density matrix from molecular orbital coefficients and occupations.
        """
        return np.dot(self.mo_coeff * self.occ, self.mo_coeff.T)

    def calculate_hirshfeld_weights(self, dm):
        """
        Calculate Hirshfeld weights for each atom using pro-molecular densities.
        """
        weights = np.zeros((self.mol.natm, self.mol.nao_nr()))
        for i, atom_density in enumerate(self.promolecular_densities):
            self.mol.set_rinv_origin(self.mol.atom_coord(i))
            atomic_density_matrix = self.mol.intor("int1e_rinv") * atom_density
            weights[i] = np.sum(atomic_density_matrix * dm, axis=1)
        total_density = np.sum(weights, axis=0)
        return weights / (total_density + 1e-12)  # Avoid division by zero

    def calculate_multipole_moments(self, weights, max_l=2):
        """
        Compute multipole moments up to order `max_l` for each atom.
        Parameters:
            weights: Hirshfeld weights for each atom.
            max_l: Maximum order of multipole moments (0=monopole, 1=dipole, 2=quadrupole, etc.).
        Returns:
            multipole_moments: Dictionary of multipole moments for each atom.
        """
        multipole_moments = {l: [] for l in range(max_l + 1)}
        for i, weight in enumerate(weights):
            moments = []
            for l in range(max_l + 1):
                self.mol.set_rinv_origin(self.mol.atom_coord(i))
                moment_matrix = self.mol.intor(f"int1e_r{l}") * weight
                moments.append(np.sum(moment_matrix))
            for l, moment in enumerate(moments):
                multipole_moments[l].append(moment)
        return multipole_moments

    def calculate_hirshfeld_properties(self, max_iterations=100, tol=1e-6, max_l=2):
        """
        Compute iterative Hirshfeld charges and multipole moments.
        Parameters:
            max_iterations: Maximum number of iterations.
            tol: Convergence tolerance for charges.
            max_l: Maximum order of multipole moments.
        Returns:
            charges: Hirshfeld charges for each atom.
            multipole_moments: Multipole moments up to order `max_l`.
        """
        n_atoms = self.mol.natm
        charges = np.zeros(n_atoms)
        prev_charges = np.ones(n_atoms) * np.inf
        iteration = 0

        dm = self.calculate_density_matrix()

        while iteration < max_iterations and np.linalg.norm(charges - prev_charges) > tol:
            iteration += 1
            prev_charges = charges.copy()

            # Update weights and calculate charges
            weights = self.calculate_hirshfeld_weights(dm)
            charges = np.einsum("ij,ij->i", weights, dm) - self.mol.atom_charges()

        if iteration >= max_iterations:
            print("Warning: Iterative Hirshfeld charges did not converge.")
        else:
            print(f"Hirshfeld charges converged in {iteration} iterations.")

        # Calculate multipole moments
        multipole_moments = self.calculate_multipole_moments(weights, max_l=max_l)
        return charges, multipole_moments

    def run(self, max_l=2):
        """
        Run the full iterative Hirshfeld calculation.
        Parameters:
            max_l: Maximum order of multipole moments.
        Returns:
            charges: Converged Hirshfeld charges.
            multipole_moments: Multipole moments up to order `max_l`.
        """
        self.load_molden_data()
        self.load_promolecular_densities()
        return self.calculate_hirshfeld_properties(max_l=max_l)

# Example usage
if __name__ == "__main__":
    molden_file = "path_to_your_molden_file.molden"  # Replace with your file path
    rad_file = "path_to_your_rad_file.rad"          # Replace with your .rad file path

    hirshfeld_calculator = IterativeHirshfeld(molden_file, rad_file)
    charges, multipole_moments = hirshfeld_calculator.run(max_l=2)

    print("Iterative Hirshfeld Charges:")
    for i, charge in enumerate(charges):
        print(f"Atom {i + 1} ({hirshfeld_calculator.mol.atom_symbol(i)}): {charge:.4f}")

    print("Multipole Moments:")
    for l, moments in multipole_moments.items():
        print(f"Order {l} Moments:")
        for i, moment in enumerate(moments):
            print(f"  Atom {i + 1}: {moment:.4f}")

