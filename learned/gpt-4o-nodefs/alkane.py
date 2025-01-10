"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is a saturated hydrocarbon with only single carbon-carbon bonds, 
    and no cyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for non-carbon and non-hydrogen elements
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):  # 1 for hydrogen, 6 for carbon
            return False, "Contains elements other than carbon and hydrogen"

    # Check for presence of double or triple bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE,):  # Only single bonds
            return False, "Contains non-single (double/triple) bonds"

    # Check for cyclic structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic structures"

    return True, "Molecule is an alkane with all single bonds and no cyclic structures"