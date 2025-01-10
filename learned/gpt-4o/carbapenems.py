"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem has a beta-lactam ring fused to a sulfur-containing five-membered ring,
    with specific stereochemistry and various substitutions.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) True if a carbapenem with reason, False with reason otherwise
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Beta-lactam pattern: 4-membered ring with nitrogen (N) and carbonyl (C=O)
    beta_lactam_pattern = Chem.MolFromSmarts("N1C(=O)C(*)C1")  # '*'-wildcard for carbon/hydrogen side chains
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"
    
    # Fused thiazoline-like sulfur-containing ring
    sulfur_ring_pattern = Chem.MolFromSmarts("C1SC2N=C1C=O")  # Simplistic approach to identify a thiazolidine structure
    if not mol.HasSubstructMatch(sulfur_ring_pattern):
        return False, "No sulfur-containing ring fused to beta-lactam found"

    # Return True if both patterns match, suggesting a carbapenem structure
    return True, "Contains the characteristic structure of a carbapenem antibiotic"

# This function can be tested on known carbapenem SMILES entries to validate its accuracy.