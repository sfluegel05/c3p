"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is characterized by a nitrogen atom bonded to two carbon atoms and one hydrogen atom,
    excluding configurations where nitrogen is part of an interfering group like nitroso, urea, or amide.

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

    # Secondary amine detection pattern with additional specificity
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1][C;$([CDATA1,#6X4]),!$(C=#[!#6])][C,$([CDATA1,#6X4]),!$(C=#[!#6])]")
        
    # Interfering group patterns
    nitroso_pattern = Chem.MolFromSmarts("[NX2]=O")
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    urea_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")
    sulfonamide_pattern = Chem.MolFromSmarts("[NX3][#16](=[OX1])=[OX1]")

    if mol.HasSubstructMatch(secondary_amine_pattern):
        # Check for interfering groups
        if any(mol.HasSubstructMatch(pattern) for pattern in [nitroso_pattern, amide_pattern, urea_pattern, sulfonamide_pattern]):
            return False, "Contains an interfering group (e.g., nitroso, urea, sulfonamide, or amide), not a secondary amine"
        
        return True, "Contains a nitrogen atom bonded to two carbon atoms and one hydrogen atom, characteristic of a secondary amine"
    
    return False, "Does not satisfy the structural requirements for a secondary amine"