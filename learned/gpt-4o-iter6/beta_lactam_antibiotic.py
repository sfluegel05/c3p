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
        Chem.MolFromSmarts("C1(C=O)N1"),  # General 4-membered ring with nitrogen
        
        # Penicillin core: thiazolidine fused with beta-lactam
        Chem.MolFromSmarts("C1([C@H])SC[C@](N1)C(=O)"),  # Thiazolidine ring with chirality
        
        # Cephalosporin core: dihydrothiazine fused with beta-lactam
        Chem.MolFromSmarts("C1([C@H])SC2=C(CN1)C(=O)"),  # Dihydrothiazine ring with chirality
        
        # Carbapenem core: unsaturated beta-lactam with sulfur
        Chem.MolFromSmarts("C1C(=C)N(C1=O)"),  # Unsaturated nitrogen-containing 4-membered ring
    ]
    
    # Pattern constraints to refine false positives
    addtional_constraints = [
        Chem.MolFromSmarts("S1C([C@@H](N1)C(=O)O")  # Constraint capturing sulfur atom often found in antibiotics's cores
    ]
    
    # Check if any of the beta-lactam patterns match
    matches_betalactam = any(mol.HasSubstructMatch(pattern) for pattern in beta_lactam_patterns)
    matches_constraints = any(not mol.HasSubstructMatch(constraint) for constraint in addtional_constraints)
    
    if matches_betalactam and not matches_constraints:
        return True, "Contains a beta-lactam ring in a known antibiotic scaffold"
    
    return False, "No beta-lactam ring found or not in a typical antibiotic context"