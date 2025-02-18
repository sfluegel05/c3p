"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene (acyclic hydrocarbon with exactly one C=C bond and formula CnH2n)
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is defined as an acyclic branched or unbranched hydrocarbon having exactly 
    one carbon-carbon double bond and the general formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an alkene, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the hydrogen count is accurate
    mol = Chem.AddHs(mol)
    
    # Check that the molecule is a strict hydrocarbon (only Carbon and Hydrogen)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than carbon and hydrogen"
    
    # Ensure the molecule is acyclic (has no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected an acyclic structure"
    
    # Count the number of carbon-carbon double bonds in the molecule
    double_bond_count = 0
    for bond in mol.GetBonds():
        # Check if the bond is a double bond between two carbons
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                double_bond_count += 1
    if double_bond_count != 1:
        return False, f"Molecule has {double_bond_count} carbon-carbon double bond(s); expected exactly 1"
    
    # Count the number of carbons and hydrogens in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    # Check if the molecular formula follows CnH2n
    if h_count != 2 * c_count:
        return False, f"Molecular formula is C{c_count}H{h_count}, which does not match CnH2n"
    
    return True, f"Molecule is a valid acyclic alkene with 1 C=C bond and formula C{c_count}H{h_count}"