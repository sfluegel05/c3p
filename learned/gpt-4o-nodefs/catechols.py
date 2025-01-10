"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is characterized by the presence of a 1,2-dihydroxybenzene moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for catechol: 1,2-dihydroxybenzene
    catechol_pattern = Chem.MolFromSmarts("Oc1ccc(O)c1")
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains catechol (1,2-dihydroxybenzene) moiety"

    return False, "Does not contain catechol moiety"

# Example usage
smiles_example = "Oc1ccccc1O"  # Catechol example
result, reason = is_catechols(smiles_example)
print(result, reason)