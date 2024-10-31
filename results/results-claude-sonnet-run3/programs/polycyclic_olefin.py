from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors


def is_polycyclic_olefin(smiles: str):
    """
    Determines if a molecule is a polycyclic olefin (polycyclic hydrocarbon with double bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic olefin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check if molecule contains only C and H atoms
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if not all(symbol in ['C', 'H'] for symbol in atoms):
        return False, "Contains non-hydrocarbon atoms"

    # Count number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 2:
        return False, "Not polycyclic (less than 2 rings)"

    # Check for presence of double bonds
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds += 1

    if double_bonds == 0:
        return False, "No double bonds present"

    return True, f"Polycyclic hydrocarbon with {num_rings} rings and {double_bonds} double bonds"
# Pr=1.0
# Recall=0.9058823529411765