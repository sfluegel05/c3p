"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic hydrocarbon with exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains only C and H atoms
    atoms = mol.GetAtoms()
    for atom in atoms:
        if atom.GetAtomicNum() not in [1, 6]:  # Only H and C allowed
            return False, "Contains non-hydrocarbon atoms"

    # Check for cycles
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings (must be acyclic)"

    # Count carbon-carbon double bonds
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                double_bonds += 1
            else:
                return False, "Contains non-carbon-carbon double bond"

    if double_bonds != 1:
        return False, f"Contains {double_bonds} carbon-carbon double bonds (must have exactly 1)"

    # Count carbons and hydrogens
    c_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 1)
    
    # Verify formula CnH2n
    expected_h = 2 * c_count
    if h_count != expected_h:
        return False, f"Formula is C{c_count}H{h_count}, expected C{c_count}H{expected_h} for alkene"

    # Additional check for total bonds
    total_bonds = sum(bond.GetBondTypeAsDouble() for bond in mol.GetBonds())
    expected_bonds = (2 * c_count - 1)  # For alkene: single bonds + one double bond
    if total_bonds != expected_bonds:
        return False, "Incorrect number of total bonds for alkene"

    return True, f"Acyclic hydrocarbon with exactly one carbon-carbon double bond (C{c_count}H{h_count})"