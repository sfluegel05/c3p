"""
Classifies: CHEBI:18310 alkane
"""
"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic, fully saturated hydrocarbon with the general formula CnH2n+2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings (must be acyclic)"

    # Check for non-carbon and non-hydrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:  # Only allow H (1) and C (6)
            return False, "Molecule contains atoms other than carbon and hydrogen"

    # Check for unsaturation (double or triple bonds)
    for bond in mol.GetBonds():
        if bond.GetBondType() not in [Chem.BondType.SINGLE]:
            return False, "Molecule contains double or triple bonds (must be fully saturated)"

    # Verify the general formula CnH2n+2
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    expected_h_count = 2 * c_count + 2
    if h_count != expected_h_count:
        return False, f"Hydrogen count ({h_count}) does not match the general formula CnH2n+2 (expected {expected_h_count})"

    return True, "Molecule is an acyclic, fully saturated hydrocarbon with the general formula CnH2n+2"