"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:38199 beta-lactam antibiotic
An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
"""
from rdkit import Chem

def is_beta_lactam_antibiotic(smiles: str):
    """
    Determines if a molecule is a beta-lactam antibiotic based on its SMILES string.

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

    # Check for beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("[C&R]2(=[O])[NR1][C&R]([C&R]2=[O])=[O]")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Additional checks for common features in beta-lactam antibiotics
    thiazolidine_pattern = Chem.MolFromSmarts("C1NCCCS1")
    dihydrothiazine_pattern = Chem.MolFromSmarts("C1NCCCS1=O")
    if mol.HasSubstructMatch(thiazolidine_pattern) or mol.HasSubstructMatch(dihydrothiazine_pattern):
        return True, "Contains beta-lactam ring and thiazolidine/dihydrothiazine rings (common in beta-lactam antibiotics)"

    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return True, "Contains beta-lactam ring and carboxylic acid group (common in beta-lactam antibiotics)"

    return True, "Contains beta-lactam ring"