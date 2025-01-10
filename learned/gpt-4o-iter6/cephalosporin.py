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
    
    # Proper beta-lactam ring pattern: 4-membered lactam with carbon-nitrogen amidic bond
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)NC1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Six-membered dihydrothiazine ring pattern: helps make it unique from penicillins
    dihydrothiazine_pattern = Chem.MolFromSmarts("C1(C(=O)N2C(S1)C2)")
    if not mol.HasSubstructMatch(dihydrothiazine_pattern):
        return False, "No dihydrothiazine ring found"

    # Ensure it does not have a penicillin-like five-membered thiazolidine ring
    penicillin_pattern = Chem.MolFromSmarts("C1CN2C(S1)C(=O)N(C2)")
    if mol.HasSubstructMatch(penicillin_pattern):
        return False, "Molecule is more penicillin-like"
    
    return True, "Molecule contains a beta-lactam and a 6-membered dihydrothiazine ring consistent with cephalosporins"

__metadata__ = {
    'chemical_class': {
        'name': 'cephalosporin',
        'definition': 'A class of beta-lactam antibiotics differing from the penicillins in having a 6-membered, rather than a 5-membered, side ring.'
    }
}