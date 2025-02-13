"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is characterized by the ester group (-C(=O)O-) where the acid part 
    is specifically octanoic acid (8 carbon chain).

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
    
    # Octanoate ester pattern: a carbonyl (C=O) with a single-bonded oxygen connected to an 8-carbon chain
    octanoate_pattern = Chem.MolFromSmarts("C(=O)O[CH2]CCCCCCC")
    if not mol.HasSubstructMatch(octanoate_pattern):
        return False, "No octanoate ester group found"
    
    return True, "Contains octanoate ester group"