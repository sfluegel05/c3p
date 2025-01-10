"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    Isoflavones are characterized by a 3-phenylchromen-4-one skeleton, often with various methoxy, hydroxyl, and glycoside substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated isoflavone pattern to account for common variations
    isoflavone_pattern = Chem.MolFromSmarts('O=c1cc(o)cc2c1ccc(c2)[cH1]c3ccccc3')
    
    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Contains isoflavone core structure (3-phenylchromen-4-one) with potential substituents"
    else:
        return False, "Does not contain isoflavone core structure"