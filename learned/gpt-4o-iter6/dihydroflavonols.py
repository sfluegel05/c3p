"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is defined as hydroxyflavanone with a hydroxyl group at position 3 of the ring structure.

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

    # Define a broader SMARTS pattern for dihydroflavonol
    # The pattern should capture the basic flavan-3-ol structure with flexibility for substituents and stereochemistry
    base_pattern = "[OH][C@]1([C@@H]([OH])[C@H](=O)C2=CC=C(C=C2)[OH])OC3=CC=CC=C31"
    dihydroflavonol_pattern = Chem.MolFromSmarts(base_pattern)
    
    # Check for basic dihydroflavonol structure
    if mol.HasSubstructMatch(dihydroflavonol_pattern):
        return True, "Contains dihydroflavonol core structure"
    
    # Test with generalized patterns accounting for variations
    generic_pattern = Chem.MolFromSmarts("O[C@H]([C@@H]([OH])C(=O)c1ccccc1)C2(OCc3ccccc32)")
    if mol.HasSubstructMatch(generic_pattern):
        return True, "Contains generalized dihydroflavonol structure with flexible stereochemistry"
    
    return False, "Does not match any dihydroflavonol backbone patterns"

# Example usage: is_dihydroflavonols("OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1")