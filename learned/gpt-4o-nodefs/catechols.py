"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol or related moiety based on its SMILES string.
    Catechol is characterized by the presence of a 1,2-dihydroxybenzene moiety directly,
    or in a form that suggests such a phenolic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a catechol or related moiety, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for catechol and related structures
    patterns = [
        Chem.MolFromSmarts("Oc1ccc(O)c1"),  # Core catechol (1,2-dihydroxybenzene)
        Chem.MolFromSmarts("Oc1c(O)cccc1"),  # To ensure overlooking by substituents
        Chem.MolFromSmarts("c1c(O)cccc1O")  # Other orientation for detecting phenolic OH
    ]

    # Check for presence of catechol or related patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains catechol or related moiety (detected)"

    return False, "Does not contain catechol or related moiety"

# Example usage
smiles_example = "Oc1ccc(O)c1"  # Catechol example
result, reason = is_catechols(smiles_example)
print(result, reason)