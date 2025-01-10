"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound contains at least one carbon-iodine (C-I) bond.

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

    # Check for presence of C-I bond
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # Check if one atom is carbon (atomic number 6) and the other is iodine (atomic number 53)
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 53) or (atom1.GetAtomicNum() == 53 and atom2.GetAtomicNum() == 6):
            return True, "Contains a carbon-iodine bond"

    return False, "No carbon-iodine bond found"

# Example usage (uncomment to test)
# example_smiles = "COCC(=O)Nc1c(I)c(C(=O)NCC(O)CO)c(I)c(C(=O)N(C)CC(O)CO)c1I"
# print(is_organoiodine_compound(example_smiles))