"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carbon-iodine bond
    ci_bond_found = False
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        
        # Check if the bond is between carbon (atomic number 6) and iodine (atomic number 53)
        if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 53) or \
           (begin_atom.GetAtomicNum() == 53 and end_atom.GetAtomicNum() == 6):
            ci_bond_found = True
            break
    
    if ci_bond_found:
        return True, "Contains at least one carbon-iodine bond"
    else:
        return False, "No carbon-iodine bond found"