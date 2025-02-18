"""
Classifies: CHEBI:18310 alkane
"""
"""
Classifies: CHEBI:18308 alkane
"""
from rdkit import Chem

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is a saturated acyclic hydrocarbon with formula CnH2n+2.

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
    
    # Add hydrogens to account for all hydrogen atoms
    mol = Chem.AddHs(mol)
    
    # Check for disconnected fragments
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) > 1:
        return False, "Multiple disconnected fragments"
    
    # Check all atoms are carbon or hydrogen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6]:
            return False, "Contains non-carbon/hydrogen atoms"
    
    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings"
    
    # Check all bonds are single bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            return False, "Contains double/triple bonds"
    
    # Count carbon and hydrogen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    # Verify molecular formula CnH2n+2
    if h_count != 2 * c_count + 2:
        return False, f"H count ({h_count}) does not match alkane formula (2*{c_count}+2 = {2*c_count+2})"
    
    return True, "Contains only carbon and hydrogen with single bonds, no cycles, and follows formula CnH2n+2"