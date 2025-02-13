"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol contains an o-diphenol component (two hydroxyl groups on adjacent carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for catechol: A benzene ring with ortho (adjacent) hydroxyl groups
    catechol_patterns = [
        Chem.MolFromSmarts("c1c(O)c(O)ccc1"),  # Basic catechol pattern
        Chem.MolFromSmarts("c1c(O)ccc(O)c1"),  # Larger aromatic system
        Chem.MolFromSmarts("c1cc(O)c(O)c2ccccc12"),  # Extended to capture more connectivity
    ]

    # Check each pattern
    for pattern in catechol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains o-diphenol component in an aromatic ring system"
    
    return False, "Does not contain o-diphenol component in an aromatic ring system"