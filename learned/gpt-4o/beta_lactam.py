"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam should have a four-membered ring containing the amide 
    bond, which involves the nitrogen and the carbonyl carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Refined SMARTS pattern for beta-lactam
    # This SMARTS pattern: "[N;H1;R]=,:[C;R](=[O;R])[C;R][C;R]" describes
    # N single-bonded within a ring -> C single-bonded -> C double-bonded to O 
    # all within a four-membered ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)N[*]1")
    
    # Check for the refined beta-lactam pattern
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a 4-membered beta-lactam ring"
    
    return False, "Does not contain a 4-membered beta-lactam ring"