"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl (CHEBI:XXXXX)
"""
from rdkit import Chem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    A para-terphenyl has a 1,4-diphenylbenzene skeleton, which is three benzene rings connected in a linear para fashion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the para-terphenyl core
    # Core is a central benzene with two benzene rings attached at para positions (1 and 4)
    # The pattern allows any substitution on the benzene rings
    para_terphenyl_core = Chem.MolFromSmarts('c1ccc(-c2ccccc2)cc1-c3ccccc3')
    
    # Check if the core structure is present
    if mol.HasSubstructMatch(para_terphenyl_core):
        return True, "Contains 1,4-diphenylbenzene skeleton"
    
    # Alternative SMARTS pattern to account for possible substitutions on the benzene rings
    # This pattern uses 'c' (aromatic) and allows any substitution on the attached benzene rings
    para_terphenyl_flex = Chem.MolFromSmarts('c1(-c2:c:c:c:c:c:2)ccc(-c3:c:c:c:c:c:3)cc1')
    if mol.HasSubstructMatch(para_terphenyl_flex):
        return True, "Contains 1,4-diphenylbenzene skeleton with possible substitutions"
    
    # If neither pattern matches
    return False, "No para-terphenyl core structure found"