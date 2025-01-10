"""
Classifies: CHEBI:22307 aldoxime
"""
from rdkit import Chem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime typically contains a carbon-nitrogen double bond (C=N) 
    and a hydroxyl group directly attached to the nitrogen (C=N-O structure).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldoxime, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible SMARTS pattern for an aldoxime, allowing for isomerism and substitution
    aldoxime_pattern = Chem.MolFromSmarts("[#6]=[N]O")  # Highlighting that N-O follows =C

    if mol.HasSubstructMatch(aldoxime_pattern):
        return True, "Contains the characteristic C=N-O structure of aldoximes"
    
    return False, "Does not contain the characteristic C=N-O structure of aldoximes"