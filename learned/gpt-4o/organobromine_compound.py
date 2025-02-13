"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound has at least one carbon-bromine (C–Br) bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all bonds in the molecule
    for bond in mol.GetBonds():
        # Get both atoms of the bond
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        
        # Check if a C-Br bond is present
        if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 35) or 
            (atom1.GetAtomicNum() == 35 and atom2.GetAtomicNum() == 6)):
            return True, "Contains at least one carbon-bromine (C–Br) bond"

    # Return failure if no C-Br bonds are found
    return False, "No carbon-bromine (C–Br) bonds found"