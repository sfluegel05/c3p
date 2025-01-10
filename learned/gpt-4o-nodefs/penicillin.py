"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillins are characterized by a beta-lactam ring and a thiazolidine ring,
    typically with a carboxyl group adjacent to the beta-lactam group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General beta-lactam pattern (four-membered lactam ring without restrictive chiral centers)
    beta_lactam_pattern = Chem.MolFromSmarts("C1(=O)NC[C@@H]1")
    
    # General thiazolidine ring pattern (five-membered sulfur-containing ring)
    thiazolidine_pattern = Chem.MolFromSmarts("C1SCCN1")
    
    # Carboxyl group pattern typically near beta-lactam
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for the core elements of the penicillin structure
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    if not mol.HasSubstructMatch(thiazolidine_pattern):
        return False, "No thiazolidine ring found"
    
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxyl group found"

    return True, "Contains beta-lactam, thiazolidine, and carboxylate group indicative of penicillins"