"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    Bile acid conjugates typically have a bile acid moiety conjugated with an amino acid or similar group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for steroid nucleus
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3(C(C(C4(C(CC(C4)O)C)C3)C)CC2)C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"
    
    # SMARTS pattern for glycine conjugation
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    
    # SMARTS pattern for taurine conjugation
    taurine_pattern = Chem.MolFromSmarts("NCCS(=O)(=O)O")
    
    # Check for glycine or taurine conjugation
    if mol.HasSubstructMatch(glycine_pattern) or mol.HasSubstructMatch(taurine_pattern):
        return True, "Bile acid conjugate with a recognized conjugation pattern"
    
    return False, "Does not have a recognized bile acid conjugation pattern"