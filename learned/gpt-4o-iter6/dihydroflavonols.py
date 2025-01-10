"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    Dihydroflavonols have a hydroxyflavanone structure with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for dihydroflavonol core: 2-phenyl-3-chromanol (OH group at C3, ketone at C4)
    core_pattern = Chem.MolFromSmarts("O[C@@H]1C[C@H](=O)[C@H](Oc2ccccc12)c3ccccc3")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match basic dihydroflavonol core structure"
    
    # Check for presence of hydroxyl groups at required positions (3rd on central ring and potential substitutions on rings)
    hydroxyl_pattern = Chem.MolFromSmarts("[$([OH]-c1ccccc1),$([OH]-c1ccc(O)cc1),$([OH]-c1ccc(O)c(O)c1)]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient hydroxyl groups at key positions"
    
    # If passes all checks, it is a dihydroflavonol
    return True, "Contains dihydroflavonol core structure and appropriate hydroxylation"

# Example usage for testing: is_dihydroflavonols("OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1")