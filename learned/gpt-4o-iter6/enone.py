"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha, beta-unsaturated ketone with the C=O function
    conjugated to a C=C double bond at the alpha, beta position, and R(4) not hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern to match enone structure
    # Target pattern: carbonyl (C=O) adjacent to alkene (C=C), preventing R4 being hydrogen
    enone_pattern = Chem.MolFromSmarts("C=CC(=O)[C;X4,!#1]")

    # Check if the molecule has the enone substructure
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains alpha, beta-unsaturated ketone (enone) structure"
    
    return False, "No enone structure found"