"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    Bile acid conjugates typically have a steroid nucleus and are conjugated with
    an amino acid or similar group.

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
    
    # General SMARTS pattern for bile acids with steroid nucleus
    steroid_pattern = Chem.MolFromSmarts("[#6]1:[#6]:[#6]:[#6]2:[#6]3:[#6]4:[#6]:1:[#6](=[#8]):[#6]:[#6]:2:[#6](=[#8]):[#6]3:[#6]:[#6]4")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No recognizable steroid nucleus found"

    # SMARTS patterns for known conjugations, expanded for common amino acids
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    taurine_pattern = Chem.MolFromSmarts("NCCS(=O)(=O)O")
    cysteine_pattern = Chem.MolFromSmarts("NCC(=O)C(O)C(S)=O")
    alanine_pattern = Chem.MolFromSmarts("NC(C)C(=O)O")
    serine_pattern = Chem.MolFromSmarts("NCC(O)C(=O)O")
    threonine_pattern = Chem.MolFromSmarts("NC(C(O)C)C(=O)O")
    
    # Check for any known conjugation
    if (mol.HasSubstructMatch(glycine_pattern) or 
        mol.HasSubstructMatch(taurine_pattern) or
        mol.HasSubstructMatch(cysteine_pattern) or 
        mol.HasSubstructMatch(alanine_pattern) or 
        mol.HasSubstructMatch(serine_pattern) or
        mol.HasSubstructMatch(threonine_pattern)):
        return True, "Bile acid conjugate with a recognized amino acid conjugation pattern"
    
    return False, "Does not have a recognized bile acid conjugation pattern"