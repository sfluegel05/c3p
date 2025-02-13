"""
Classifies: CHEBI:32878 alkene
"""
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

    # Check if molecule is cyclic
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings - must be acyclic"

    # Check if molecule contains only C and H atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:  # 1=H, 6=C
            return False, "Contains non-hydrocarbon atoms"

    # Count carbon-carbon double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    n_double_bonds = len(double_bond_matches)
    
    if n_double_bonds == 0:
        return False, "No carbon-carbon double bonds found"
    elif n_double_bonds > 1:
        return False, "Contains multiple double bonds"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for any other types of unsaturation (triple bonds, etc.)
    # by comparing bond counts with what's expected for an alkene
    total_bonds = sum(bond.GetBondTypeAsDouble() for bond in mol.GetBonds())
    expected_bonds = (2 * c_count - 1)  # For alkene: single bonds + one double bond
    
    if total_bonds != expected_bonds:
        return False, "Contains unexpected bonding pattern"

    return True, f"Acyclic hydrocarbon with one C=C double bond"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:32878',
        'name': 'alkene',
        'definition': 'An acyclic branched or unbranched hydrocarbon having one '
                     'carbon-carbon double bond and the general formula CnH2n. '
                     'Acyclic branched or unbranched hydrocarbons having more than '
                     'one double bond are alkadienes, alkatrienes, etc.',
        'parents': ['CHEBI:33641']
    }
}