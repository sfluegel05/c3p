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

    # More accurate beta-lactam pattern considering chiral centers
    beta_lactam_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([NX3])C(=O)[CX3])S[C@H]2N1C([C@@H]2C)=O")
    
    # Thiazolidine pattern including chiral specificity
    thiazolidine_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([N])SC2)C2C1")
    
    # Carboxylate group being considered more generally near beta-lactam
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for elements of the penicillin structure
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    if not mol.HasSubstructMatch(thiazolidine_pattern):
        return False, "No thiazolidine ring found"
    
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxyl group found"

    return True, "Contains beta-lactam, thiazolidine, and carboxylate group indicative of penicillins"