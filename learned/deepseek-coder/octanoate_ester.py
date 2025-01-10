"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:75548 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the octanoate ester pattern: C(=O)O-R where R is any group
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O[!H]")
    
    # Check if the pattern is present in the molecule
    if mol.HasSubstructMatch(octanoate_pattern):
        return True, "Contains octanoate ester group (CCCCCCCC(=O)OR)"
    else:
        return False, "No octanoate ester group found"