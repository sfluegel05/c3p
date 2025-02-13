"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is characterized by a nitrogen atom bonded to two carbon atoms and one hydrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved secondary amine SMARTS pattern
    # Detect nitrogen attached to any two carbon atoms and one hydrogen atom
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1]([#6])[#6]")
    
    if mol.HasSubstructMatch(secondary_amine_pattern):
        
        # Verify that the nitrogen is not part of a known interfering group
        # E.g., Nitroso, amide, etc.
        nitroso_pattern = Chem.MolFromSmarts("[NX2]=O")
        amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
        
        if mol.HasSubstructMatch(nitroso_pattern) or mol.HasSubstructMatch(amide_pattern):
            return False, "Contains an interfering group (e.g., nitroso or amide), not a secondary amine"
    
        return True, "Contains a nitrogen atom bonded to two carbon atoms and one hydrogen atom, characteristic of a secondary amine"
    
    return False, "Does not satisfy the structural requirements for a secondary amine"