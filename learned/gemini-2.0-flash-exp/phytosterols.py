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

    # Define a more flexible steroid core pattern
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C]([C]4[C]([C]12)CC[C]34)[C]")
    if not mol.HasSubstructMatch(steroid_core_pattern):
            return False, "Does not have the core steroid structure"
            
    # Check for hydroxyl group at position 3 (we need to be precise here)
    hydroxyl_pattern = Chem.MolFromSmarts("[C]1[C]([OH])[C]2[C]3[C]([C]4[C]([C]12)CC[C]34)[C]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Does not have hydroxyl group at position 3"
    
    # Check for a side chain at position 17
    sidechain_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C]([C]4[C]([C]12)CC[C]34)[C]-[C]")
    if not mol.HasSubstructMatch(sidechain_pattern):
            return False, "Does not have a side chain at position 17"


    # Check for double bonds (optional)
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    if mol.HasSubstructMatch(double_bond_pattern):
        return True, "Has core steroid structure with side chain, hydroxyl group at 3 position and at least one double bond"
    

    return True, "Has core steroid structure with side chain and hydroxyl group at 3 position"