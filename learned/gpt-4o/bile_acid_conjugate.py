"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    
    A bile acid conjugate is defined as any bile acid conjugated to a hydrophilic or charged functional group, 
    typically including glycine, taurine, sulfate, glucuronic acid, and others.
    The core structure should exhibit the steroid skeleton of bile acids, with acceptable conjugation.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved steroid skeleton pattern (includes stereo centers and features of bile acids)
    steroid_skeleton_pattern = Chem.MolFromSmarts("C1[C@@H]2C[C@H]3[C@H]([C@@H](CC(=O)O)C)CC[C@]4(C)C(=O)[C@H]([C@@]12C)C[C@@H]4O")
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroidal skeleton typical of bile acids found"

    # Conjugation moiety patterns (additional conjugates are checked)
    conjugate_patterns = {
        'glycine': Chem.MolFromSmarts("NCC(=O)O"),
        'taurine': Chem.MolFromSmarts("NCCS(=O)(=O)O"),
        'sulfate': Chem.MolFromSmarts("OS(=O)(=O)O"),
        'glucuronic': Chem.MolFromSmarts("OC[C@@H]1O[C@H](C(O)[C@H](O)C1=O)"),
        'alanine': Chem.MolFromSmarts("NC(C)C(=O)O"),
        'valine': Chem.MolFromSmarts("NC(C(C)C)C(=O)O"),
        'serine': Chem.MolFromSmarts("NC(CO)C(=O)O"),
        'phenylalanine': Chem.MolFromSmarts("NC(CC1=CC=CC=C1)C(=O)O"),
        'arginine': Chem.MolFromSmarts("NC(CCCNC(N)=N)C(=O)O"),
    }

    # Check for conjugation patterns on the molecule
    for conjugate_name, pattern in conjugate_patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains steroid skeleton with {conjugate_name} conjugate"
    
    # If none of the conjugate is detected
    return False, "No recognizable bile acid conjugate structure found"