"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compounds
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine based on its SMILES string.
    An organobromine compound contains at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbon-bromine bonds in the molecule
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        atomic_nums = sorted([atom1.GetAtomicNum(), atom2.GetAtomicNum()])
        if atomic_nums == [6, 35]:  # 6 is atomic number for carbon, 35 for bromine
            return True, "Contains at least one carbon-bromine bond"

    return False, "No carbon-bromine bond found"