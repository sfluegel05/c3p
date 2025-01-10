"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic contains a beta-lactam ring: a 4-membered ring with one nitrogen atom.
    We consider common antibiotic classes like penicillins, cephalosporins, carbapenems, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define broader SMARTS patterns for beta-lactam rings
    # This accounts for different stereochemistry and larger systems:
    beta_lactam_patterns = [
        # Beta-lactam core: 4-membered ring with a nitrogen and carbonyl
        Chem.MolFromSmarts("C1CNC(=O)1"),
        
        # Penicillin-like core: 5-membered thiazolidine fused with beta-lactam ring
        Chem.MolFromSmarts("C1CNC(=O)1C"),
        
        # Cephalosporin-like core: 6-membered dihydrothiazine fused with beta-lactam
        Chem.MolFromSmarts("C1CNC(=O)1C2CS2"),
    ]
    
    # Check if any of the beta-lactam patterns match
    for pattern in beta_lactam_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a beta-lactam ring in a common antibiotic scaffold"
    
    return False, "No beta-lactam ring found or not in a typical antibiotic context"