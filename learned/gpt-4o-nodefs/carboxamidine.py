"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group based on its SMILES string.
    A carboxamidine is characterized by the structure R-N=C(NR)NR, where R can be a hydrogen 
    or any other substituent. The function attempts to recognize various possible carboxamidine
    structures present within the molecule.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a carboxamidine group, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an array of potential SMARTS patterns for carboxamidine
    # Adding more specific patterns based on examples and avoiding matches with certain structures
    carboxamidine_patterns = [
        Chem.MolFromSmarts("N=C(N)N"),  # Basic carboxamidine form
        Chem.MolFromSmarts("N=C(N)N"),  # Simple form allowing symmetry
        Chem.MolFromSmarts("[NX3][CX3]=[NX2][NX3]"),  # Allowing for varied substitutions
        # Add more detailed patterns derived from examples and known structures
    ]

    # Check each pattern
    for pattern in carboxamidine_patterns:
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carboxamidine group matching the structure"

    return False, "Does not contain recognizable carboxamidine groups"