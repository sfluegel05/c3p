"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha, beta-unsaturated ketone where the C=O function is conjugated to a C=C double bond
    at the alpha, beta position, and R(4) is not hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecular object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for an enone: C=C-C(=O)-C (R4 shouldn't be hydrogen)
    enone_pattern = Chem.MolFromSmarts("[CX3](=O)[CX3;!R]=[CX3][CX3;!Ha]") 

    # Check for the enone pattern in the molecule
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains alpha, beta-unsaturated ketone (enone) structure"

    return False, "No enone structure found"