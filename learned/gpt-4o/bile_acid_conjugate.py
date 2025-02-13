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

    # Revised steroid skeleton pattern
    steroid_skeleton_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC(C4)C3CCC12")  # Generalized cyclopenta[a]phenanthrene pattern
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroidal skeleton typical of bile acids found"

    # Conjugation moiety patterns
    conjugate_patterns = {
        'glycine': Chem.MolFromSmarts("NCC(=O)O"),
        'taurine': Chem.MolFromSmarts("NCCS(=O)(=O)O"),
        'sulfate': Chem.MolFromSmarts("O[S](=O)(=O)O"),
        'glucuronic': Chem.MolFromSmarts("OC[C@@H]1O[C@H](C(O)[C@H](O)C1=O)"),
        'alanine': Chem.MolFromSmarts("NC(C)C(=O)O"),
        'valine': Chem.MolFromSmarts("NC(C(C)C)C(=O)O"),
        'serine': Chem.MolFromSmarts("NC(CO)C(=O)O"),
        'phenylalanine': Chem.MolFromSmarts("NC(CC1=CC=CC=C1)C(=O)O"),
        'arginine': Chem.MolFromSmarts("NC(CCCNC(N)=N)C(=O)O"),
        # Add more conjugation patterns as needed
    }

    # Check for conjugation
    for conjugate_name, pattern in conjugate_patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains steroid skeleton with {conjugate_name} conjugate"

    # If it reaches here, no conjugation was detected
    return False, "No recognizable bile acid conjugate structure found"