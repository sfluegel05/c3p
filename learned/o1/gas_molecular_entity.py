"""
Classifies: CHEBI:138675 gas molecular entity
"""
# gas_molecular_entity.py
"""
Classifies: CHEBI gas molecular entity
"""
from rdkit import Chem

def is_gas_molecular_entity(smiles: str):
    """
    Determines if a molecule is a main group molecular entity that is gaseous at STP.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a gas molecular entity, False otherwise.
        str: Reason for classification.
    """
    from rdkit.Chem import rdMolDescriptors

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for net charge; exclude ions
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule is an ion with non-zero net charge, unlikely to be gaseous at STP"

    # List of atomic numbers for noble gases
    noble_gases = [2, 10, 18, 36, 54, 86]  # He, Ne, Ar, Kr, Xe, Rn
    # List of diatomic gases
    diatomic_gases_smiles = ['[H][H]', 'N#N', 'O=O', 'F[F]', 'ClCl']
    # List of common gaseous inorganic molecules
    common_gaseous_inorganics_smiles = [
        'O=C=O',     # Carbon dioxide
        '[C-]#[O+]', # Carbon monoxide
        'N#N',       # Nitrogen
        'O=O',       # Oxygen
        'Cl[H]',     # Hydrogen chloride
        'I[H]',      # Hydrogen iodide
        '[O][O]',    # Triplet oxygen
        'Br[H]',     # Hydrogen bromide
        '[H]N([H])[H]', # Ammonia
        '[H]F',      # Hydrogen fluoride
        '[O-][O+]=O',# Ozone
    ]

    # Check if the molecule is a noble gas atom
    if mol.GetNumAtoms() == 1:
        atomic_num = mol.GetAtomWithIdx(0).GetAtomicNum()
        if atomic_num in noble_gases:
            return True, f"Molecule is a noble gas atom ({mol.GetAtomWithIdx(0).GetSymbol()}), gaseous at STP"

    # Check if the molecule matches known diatomic gases or common gaseous inorganics
    all_gaseous_smiles = diatomic_gases_smiles + common_gaseous_inorganics_smiles
    for gas_smiles in all_gaseous_smiles:
        gas_mol = Chem.MolFromSmiles(gas_smiles)
        if mol.GetNumAtoms() == gas_mol.GetNumAtoms() and Chem.MolToSmiles(mol, isomericSmiles=True) == Chem.MolToSmiles(gas_mol, isomericSmiles=True):
            return True, "Molecule is a known gaseous inorganic compound at STP"

    # Implement the Joback method to estimate boiling point
    # Group contribution values for the Joback method
    joback_groups = {
        ('C', 4): {'Tb': 23.58},  # -CH3
        ('C', 3): {'Tb': 22.88},  # >CH2
        ('C', 2): {'Tb': 21.74},  # >CH-
        ('C', 1): {'Tb': 18.90},  # >C<
        ('H', 0): {'Tb': 0.00},   # Hydrogen (does not contribute)
        ('O', 2): {'Tb': 17.40},  # -OH (alcohol)
        ('O', 1): {'Tb': 9.20},   # Ether oxygen
        ('N', 3): {'Tb': 13.90},  # -NH2
        ('F', 1): {'Tb': 9.50},   # -F
        ('Cl', 1): {'Tb': 19.50}, # -Cl
        ('Br', 1): {'Tb': 26.30}, # -Br
        ('I', 1): {'Tb': 33.80},  # -I
        ('S', 2): {'Tb': 22.60},  # -SH
        # Add more groups as needed
    }

    # Function to get Joback group counts
    from collections import defaultdict

    def get_joback_group_counts(mol):
        group_counts = defaultdict(int)
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            valence = atom.GetTotalValence()
            hydrogens = atom.GetTotalNumHs()
            bonds = atom.GetDegree()
            charge = atom.GetFormalCharge()
            hybridization = atom.GetHybridization()

            key = (symbol, bonds)
            if key in joback_groups:
                group_counts[key] += 1
            else:
                # Approximate unknown groups
                key = (symbol, 0)
                group_counts[key] += 1
        return group_counts

    # Estimate boiling point using Joback method
    group_counts = get_joback_group_counts(mol)
    Tb = 198.0  # Base temperature
    for key, count in group_counts.items():
        contribution = joback_groups.get(key, {}).get('Tb', 0) * count
        Tb += contribution

    # Adjust for molecular structure (cyclic compounds)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 0:
        Tb -= 33.0 * num_rings  # Approximate adjustment

    Tb_K = Tb  # Boiling point in Kelvin
    Tb_C = Tb_K - 273.15  # Convert to Celsius

    # Classify based on estimated boiling point
    if Tb_C <= 0:
        reason = f"Estimated boiling point is {Tb_C:.1f}°C, molecule is gaseous at STP"
        return True, reason
    else:
        return False, f"Estimated boiling point is {Tb_C:.1f}°C, molecule is not gaseous at STP"