"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    A cephalosporin is characterized by having a beta-lactam ring and a 6-membered
    dihydrothiazine ring fused together with variable side chains.

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
    
    # Check for beta-lactam ring (4-membered lactam with carbon-nitrogen amidic bond)
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)N[C@@H]1")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No beta-lactam ring found"

    # Check for 6-membered dihydrothiazine ring
    # Update the SMARTS pattern to more accurately capture the ring
    dihydrothiazine_pattern = Chem.MolFromSmarts("C1(C)C(=O)N2[C@H](SC1)C2")
    if not mol.HasSubstructMatch(dihydrothiazine_pattern):
        return False, "No dihydrothiazine ring found or incorrect stereochemistry"

    # Successful match for cephalosporin scaffold
    return True, "Molecule contains a beta-lactam and a 6-membered dihydrothiazine ring consistent with cephalosporins"

__metadata__ = {
    'chemical_class': {
        'name': 'cephalosporin',
        'definition': 'A class of beta-lactam antibiotics differing from the penicillins in having a 6-membered, rather than a 5-membered, side ring.'
    }
}