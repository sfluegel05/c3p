"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.
    A beta-lactam antibiotic typically contains a beta-lactam ring: a 4-membered ring with one nitrogen atom.
    We also consider scaffolds characteristic of known antibiotics, such as penicillins, cephalosporins, and carbapenems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam antibiotic, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for beta-lactam antibiotics
    beta_lactam_patterns = [
        Chem.MolFromSmarts("C1CNC(=O)1"),  # General beta-lactam ring
        Chem.MolFromSmarts("C1[C@@H](N1)C(=O)[C@H]"),  # Stereochemically specific beta-lactam
        
        # Penicillin core: thiazolidine fused with beta-lactam
        Chem.MolFromSmarts("C1([C@H])SC[C@](N1)C(=O)"),  # Thiazolidine ring with chirality
        
        # Cephalosporin core: dihydrothiazine fused with beta-lactam
        Chem.MolFromSmarts("C1([C@H])SC2CCN1C2=O"),  # Dihydrothiazine ring with chirality
        
        # Carbapenem core: unsaturated beta-lactam with sulfur
        Chem.MolFromSmarts("C1C(=C)N(C1=O)C(=O)"),
    ]
    
    # Check if any of the beta-lactam patterns match
    for pattern in beta_lactam_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a beta-lactam ring in a known antibiotic scaffold"
    
    return False, "No beta-lactam ring found or not in a typical antibiotic context"