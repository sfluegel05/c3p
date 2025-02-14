"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound contains at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flag to track presence of C-F bonds
    contains_c_f = False

    # Traverse through bonds to check for C-F bonds
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # Check if we've found a carbon-fluorine bond
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 9) or (atom1.GetAtomicNum() == 9 and atom2.GetAtomicNum() == 6):
            contains_c_f = True
            break
    
    if contains_c_f:
        return True, "Molecule contains at least one carbon-fluorine bond"
    else:
        return False, "No carbon-fluorine bond found in the molecule"