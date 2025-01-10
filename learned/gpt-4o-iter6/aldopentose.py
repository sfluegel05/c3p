"""
Classifies: CHEBI:33916 aldopentose
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.

    An aldopentose is defined as a pentose (5-carbon sugar) with an aldehyde
    group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldopentose, False otherwise.
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if there are exactly 5 carbon atoms (pentose sugar)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbon atoms, found {c_count}"

    # Check for the presence of an aldehyde group (R-CHO)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"
    
    return True, "Contains a 5-carbon backbone with an aldehyde group at one end"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26977',
        'name': 'aldopentose',
        'definition': 'A pentose with a (potential) aldehyde group at one end.',
        'parents': ['CHEBI:4705']  # Assuming a parent class in a real ontology
    },
    'config': {
        'llm_model_name': 'your-model',
        'f1_threshold': 0.8,
        'max_attempts': 5
    },
    'message': None,
    'attempt': 0,
    'success': True
}