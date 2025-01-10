"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is characterized by a nitrogen atom bonded to two carbon atoms and one hydrogen atom,
    excluding configurations where nitrogen is part of nitroso, amide, urea, or sulfonamide groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the secondary amine SMARTS pattern
    # [NX3] specifies a 3-connected (secondary amine) nitrogen
    # [C] specifies any carbon atom
    # Adjust pattern to detect nitrogen bonded to exactly two carbons and one other atom (possibly H)

    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;R0]([C])[C]")

    # Patterns for interfering groups (exclude these)
    nitroso_pattern = Chem.MolFromSmarts("[#7][#16](=[OX1])=[OX1]")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    urea_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")
    sulfonamide_pattern = Chem.MolFromSmarts("[NX3][#16](=[OX1])=[OX1]")

    # Check for secondary amine substructure match
    if mol.HasSubstructMatch(secondary_amine_pattern):
        # Check for any interfering functional groups
        if any(mol.HasSubstructMatch(pattern) for pattern in [nitroso_pattern, amide_pattern, urea_pattern, sulfonamide_pattern]):
            return False, "Contains interfering groups (e.g., nitroso, urea, amide, sulfonamide), not a secondary amine"
        
        return True, "Contains a nitrogen atom bonded to two carbon atoms, characteristic of a secondary amine"
    
    return False, "Does not satisfy the structural requirements for a secondary amine"