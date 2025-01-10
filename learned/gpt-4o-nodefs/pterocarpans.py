"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans have a core structure of 6a,11a-dihydrobenzofuro[3,2-c]chromene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool, str: True and reason if molecule is a pterocarpan, False and reason otherwise
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the pterocarpan core structure
    pterocarpan_pattern = Chem.MolFromSmarts("C12=C3C(CC4=CC=CC=C4O3)OC56=C2C=CC(O)=C(C1=O)C56")
    
    # Check if molecule contains the pterocarpan core structure
    if mol.HasSubstructMatch(pterocarpan_pattern):
        return True, "Contains the pterocarpan core structure"
    
    return False, "Does not contain the pterocarpan core structure"