"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid skeleton pattern (generalized cyclopenta[a]phenanthrene pattern)
    steroid_skeleton_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4C3(CCC4)C")
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroidal skeleton typical of bile acids found"

    # Conjugation moiety patterns
    # Glycine pattern
    glycine_pattern = Chem.MolFromSmarts("NC(C)C(=O)O")
    # Taurine pattern
    taurine_pattern = Chem.MolFromSmarts("NCCS(=O)(=O)O")
    # Glucuronic acid pattern
    glucuronic_pattern = Chem.MolFromSmarts("OC[C@@H]1O[C@H](C(O)[C@H](O)C1=O)")
    # Sulfuric acid pattern (for sulfated bile acids)
    sulfate_pattern = Chem.MolFromSmarts("[O;R0][S;R0](=O)(=O)O")
    
    # Check for glycine, taurine, glucuronic acid, and sulfate conjugation
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Contains steroid skeleton with glycine conjugate"
    elif mol.HasSubstructMatch(taurine_pattern):
        return True, "Contains steroid skeleton with taurine conjugate"
    elif mol.HasSubstructMatch(glucuronic_pattern):
        return True, "Contains steroid skeleton with glucuronic acid conjugate"
    elif mol.HasSubstructMatch(sulfate_pattern):
        return True, "Contains steroid skeleton with sulfate conjugate"
    
    return False, "No recognizable bile acid conjugate structure found"

__metadata__ = {   'chemical_class': {   'id': 'custom',
                          'name': 'bile acid conjugate',
                          'definition': 'Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule.'}}