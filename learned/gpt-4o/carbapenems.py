"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    Carbapenems are characterized by a bicyclic core structure consisting of a 4-membered
    beta-lactam ring fused to a 5-membered ring containing sulfur, typically with stereochemistry
    and substitutions at specific positions.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) True for carbapenem with reason, False with reason otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a beta-lactam ring (4-membered ring with an -N= and -C=O)
    beta_lactam_pattern = Chem.MolFromSmarts("[NX3;H1]1C(=O)C[CH]1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Check for the presence of a 5-membered ring fused that contains sulfur
    sulfur_fused_pattern = Chem.MolFromSmarts("C1[CH]S[CH]2OC(=O)C[CH]21")  # basic carbapenem core
    if not mol.HasSubstructMatch(sulfur_fused_pattern):
        return False, "No fused sulfur-containing 5-membered ring found"

    # Optionally check for substituents in typical positions if needed given specific use
    
    return True, "Contains the characteristic structure of a carbapenem antibiotic"