"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: gamma-lactone (5-membered lactone ring)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule contains a gamma-lactone (5-membered lactone ring).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a gamma-lactone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for gamma-lactone:
    # 5-membered ring with ester group (-O-C(=O)-)
    # The numbers ensure the atoms form a 5-membered ring
    gamma_lactone_pattern = Chem.MolFromSmarts("[O;R1]1-[C;R1]-[C;R1]-[C;R1]-[C;R1](=[O;R1])1")
    
    # Alternative pattern that might catch other variants
    gamma_lactone_pattern2 = Chem.MolFromSmarts("[O;R1]1-[C;R1]-[C;R1]=[C;R1]-[C;R1](=[O;R1])1")
    
    # Check for matches
    matches = mol.GetSubstructMatches(gamma_lactone_pattern)
    matches2 = mol.GetSubstructMatches(gamma_lactone_pattern2)
    
    all_matches = len(matches) + len(matches2)
    
    if all_matches > 0:
        if all_matches == 1:
            return True, "Contains a gamma-lactone (5-membered lactone ring)"
        else:
            return True, f"Contains {all_matches} gamma-lactone rings"
    
    # Additional check for potential tautomers or variants
    variant_pattern = Chem.MolFromSmarts("[O;R1]1-[C;R1]=[C;R1]-[C;R1]-[C;R1](=[O;R1])1")
    variant_matches = mol.GetSubstructMatches(variant_pattern)
    
    if len(variant_matches) > 0:
        return True, "Contains a gamma-lactone (variant form)"
        
    return False, "No gamma-lactone substructure found"