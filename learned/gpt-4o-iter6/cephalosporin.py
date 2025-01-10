"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    A cephalosporin is characterized by having a beta-lactam ring and a 6-membered
    dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cephalosporin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-lactam ring pattern (4-membered cyclic amide)
    beta_lactam_pattern = Chem.MolFromSmarts("C1CN([C@@H]1)C(=O)")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Define the six-membered dihydrothiazine ring pattern
    dihydrothiazine_pattern = Chem.MolFromSmarts("C1SC([NX3H,OX2H1])C(=O)N([C@H]1)")
    if not mol.HasSubstructMatch(dihydrothiazine_pattern):
        return False, "No dihydrothiazine ring found"

    # Ensure it does not have a penicillin-like structure
    penicillin_pattern = Chem.MolFromSmarts("C1CN2C(S1)C(=O)N(C2)C")
    if mol.HasSubstructMatch(penicillin_pattern):
        return False, "Molecule is more penicillin-like"

    return True, "Molecule contains a beta-lactam and a dihydrothiazine ring consistent with cephalosporins"

__metadata__ = {
    'chemical_class': {
        'name': 'cephalosporin',
        'definition': 'A class of beta-lactam antibiotics differing from the penicillins in having a 6-membered, rather than a 5-membered, side ring.'
    }
}