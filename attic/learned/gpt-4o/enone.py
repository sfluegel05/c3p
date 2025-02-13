"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha, beta-unsaturated ketone of the form R(1)R(2)C=CR(3)-C(=O)R(4),
    where R(4) is not a hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the enone pattern
    # enone SMARTS pattern: C=C-C=O with an assurance of R(4) not being hydrogen
    enone_pattern = Chem.MolFromSmarts("C=CC(=O)[C;!#1]")
    
    # Check if the molecule contains the enone substructure
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains alpha, beta-unsaturated ketone structure"

    return False, "Does not contain the enone pattern"