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

    # Define the beta-lactam ring pattern: a four-membered ring with a nitrogen and a carbonyl
    beta_lactam_pattern = Chem.MolFromSmarts("[N;R1][C;R1][C;R1][C;R1]=O")
    
    # Check if the pattern matches
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a four-membered ring with an amide bond (beta-lactam ring)"
    else:
        return False, "No four-membered ring with an amide bond found"