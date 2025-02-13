"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxyl group at position 3 of the C-ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for flavanone backbone with 3-hydroxy group
    flavanone_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)c2ccc(O)cc2:C1=O")
    if mol.HasSubstructMatch(flavanone_pattern):
        return True, "Contains flavanone backbone with 3-hydroxy group"
    
    # Check alternative orientation for stereo chemistry
    flavanone_pattern_alternative = Chem.MolFromSmarts("O[C@@H]1[C@@H](O)c2ccc(O)cc2:C1=O")
    if mol.HasSubstructMatch(flavanone_pattern_alternative):
        return True, "Contains flavanone backbone with 3-hydroxy group (alternative stereochemistry)"
    
    return False, "Does not contain the 3-hydroxyflavanone structure"

# Example usage: is_dihydroflavonols("OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1")