"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: CHEBI:37577 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is any ester of carbamic acid or its N-substituted derivatives,
    characterized by the functional group -O-C(=O)-N-.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for carbamate ester functional group
    carbamate_smarts = "[OX2H0][CX3](=O)[NX3,NX4]"
    carbamate_pattern = Chem.MolFromSmarts(carbamate_smarts)
    if carbamate_pattern is None:
        return (None, None)

    # Find substructure matches
    matches = mol.GetSubstructMatches(carbamate_pattern)
    if matches:
        num_matches = len(matches)
        return True, f"Contains {num_matches} carbamate ester functional group(s)"
    else:
        return False, "Does not contain carbamate ester functional group"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:37577',
        'name': 'carbamate ester',
        'definition': 'Any ester of carbamic acid or its N-substituted derivatives.',
        'parents': ['CHEBI:51200']
    }
}