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
    A carboxamidine can generally be characterized by the structure R-N=C(NR)NR,
    where R can be a hydrogen or any other substituent group, and may involve variable
    oxidation states or additional substituents.

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

    # Define potential SMARTS patterns for carboxamidine
    carboxamidine_patterns = [
        Chem.MolFromSmarts("N=C(N)N"),  # Classic form
        Chem.MolFromSmarts("[#7]C(=N-[!#1])N"),  # More generalized pattern with atoms placeholder
        Chem.MolFromSmarts("N=C(N)(N)"),  # Allowing for symmetric forms
    ]

    # Verify each pattern
    for pattern in carboxamidine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carboxamidine group matching one of the generalized patterns"

    return False, "Does not contain recognizable carboxamidine groups"