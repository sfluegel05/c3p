"""
Classifies: CHEBI:27933 beta-lactam antibiotic
"""
"""
Classifies: CHEBI:38199 beta-lactam antibiotic
An organonitrogen heterocyclic antibiotic that contains a beta-lactam ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Check for penicillin-specific features
    thiazolidine_pattern = Chem.MolFromSmarts("C1NC(=O)CS1")
    if mol.HasSubstructMatch(thiazolidine_pattern):
        return True, "Contains beta-lactam ring and thiazolidine ring (penicillin)"

    # Check for cephalosporin-specific features
    dihydrothiazine_pattern = Chem.MolFromSmarts("C1NC(=O)CS1=O")
    furan_pattern = Chem.MolFromSmarts("c1occc1")
    if mol.HasSubstructMatch(dihydrothiazine_pattern) or mol.HasSubstructMatch(furan_pattern):
        return True, "Contains beta-lactam ring and dihydrothiazine ring or furan ring (cephalosporin)"

    # Check for other common features
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    amino_pattern = Chem.MolFromSmarts("N")
    if mol.HasSubstructMatch(carboxyl_pattern) and mol.HasSubstructMatch(amino_pattern):
        return True, "Contains beta-lactam ring, carboxylic acid group, and amino group (common in beta-lactam antibiotics)"

    return True, "Contains beta-lactam ring (beta-lactam antibiotic)"