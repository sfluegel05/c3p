"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester contains an ester linkage with an octanoic acid part (8-carbon chain).
    
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

    # Define ester bond pattern including variable chain on either side
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for basic ester bond
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond found"

    # Ensure there is an 8-carbon chain linked to the ester group
    # The pattern assumes linear CH2 groups in the chain (R-C(=O)O-R')
    octanoate_chain_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    
    # Check for octanoate specifically attached to an oxygen
    if not mol.HasSubstructMatch(octanoate_chain_pattern):
        return False, "No octanoate chain found associated with ester bond"

    # If all checks pass, it is an octanoate ester
    return True, "Molecule contains an octanoate ester group"