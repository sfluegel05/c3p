"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:26680 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols that occur in plants and have a core steroid structure, a hydroxyl at position 3 and side chains
    with some double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Simplified steroid core pattern (tetracyclic ring system) - allows for double bonds and chirality
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C][C]3[C]([C]4[C]([C]1[C]2)CC[C]34)[C]")
    if not mol.HasSubstructMatch(steroid_core_pattern):
            return False, "Does not have the core steroid structure"

    # Check for side chain at C17, and check for a C attached to C17 of the core
    sidechain_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C][C]3[C]([C]4[C]([C]1[C]2)CC[C]34)[C]-[C]")
    if not mol.HasSubstructMatch(sidechain_pattern):
        return False, "Does not have the side chain at C17"
    
    # Check for hydroxyl at position 3
    hydroxyl_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C][C]3[C]([C]4[C]([C]1[C]2)CC[C]34)[C][OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Does not have hydroxyl at C3"

    # Optional check for double bonds
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    if mol.HasSubstructMatch(double_bond_pattern):
        return True, "Has core steroid structure with side chain, hydroxyl group and at least one double bond"
    
    return True, "Has core steroid structure with side chain and hydroxyl group"