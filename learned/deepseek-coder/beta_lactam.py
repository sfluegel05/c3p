"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam is a lactam in which the amide bond is contained within a four-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible beta-lactam ring pattern: a four-membered ring with a nitrogen and a carbonyl
    # This pattern allows for other atoms in the ring and accounts for bicyclic systems
    beta_lactam_pattern = Chem.MolFromSmarts("[N;R1]1[C;R1][C;R1][C;R1]1=O")
    
    # Check if the pattern matches
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a four-membered ring with an amide bond (beta-lactam ring)"
    
    # Define an alternative pattern for bicyclic systems where the beta-lactam ring is fused to another ring
    beta_lactam_bicyclic_pattern = Chem.MolFromSmarts("[N;R1]1[C;R1][C;R1][C;R1]1=O.[C;R1]2[C;R1][C;R1][C;R1]2")
    
    # Check if the alternative pattern matches
    if mol.HasSubstructMatch(beta_lactam_bicyclic_pattern):
        return True, "Contains a four-membered ring with an amide bond (beta-lactam ring) in a bicyclic system"
    
    # Define a pattern for beta-lactams with sulfur in the ring (e.g., penicillins)
    beta_lactam_sulfur_pattern = Chem.MolFromSmarts("[N;R1]1[C;R1][S;R1][C;R1]1=O")
    
    # Check if the sulfur-containing pattern matches
    if mol.HasSubstructMatch(beta_lactam_sulfur_pattern):
        return True, "Contains a four-membered ring with an amide bond and sulfur (beta-lactam ring)"
    
    # If no pattern matches, return False
    return False, "No four-membered ring with an amide bond found"