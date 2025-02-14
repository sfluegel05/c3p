"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has a decanoyl group (10 carbons) attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the decanoyl group attached to an ester oxygen
    decanoyl_ester_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)O")

    # Check for the decanoate ester substructure
    if mol.HasSubstructMatch(decanoyl_ester_pattern):
         return True, "Contains a decanoyl group connected via an ester bond"
    
    return False, "No decanoate ester core found"