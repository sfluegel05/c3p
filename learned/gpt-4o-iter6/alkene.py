"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic branched or unbranched hydrocarbon with one C=C double bond,
    following the formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of exclusively carbon and hydrogen
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H']:
            # Allow some non-C/H component presence as long as double bonds remain unaffected
            return False, f"Molecule contains a non-CH element: {atom.GetSymbol()}"
    
    # Check for cycles: alkenes must be acyclic in terms of double bonds
    ring_info = mol.GetRingInfo()
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.IsInRing():
            return False, "C=C bond is in a ring; alkenes must have acyclic double bonds" 
    
    # Count carbon-carbon double bonds
    double_bond_ccount = sum(1 for bond in mol.GetBonds() 
                             if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
                             bond.GetBeginAtom().GetSymbol() == 'C' and 
                             bond.GetEndAtom().GetSymbol() == 'C')
    
    # Alkenes must have exactly one acyclic C=C double bond
    if double_bond_ccount != 1:
        return False, f"Expected 1 acyclic C=C double bond, found {double_bond_ccount}"
    
    # Count carbons and hydrogens, including implicit hydrogens
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    num_hydrogens = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    
    # Verify the formula CnH2n
    if num_hydrogens != 2 * num_carbons:
        return False, f"Formula check failed: CnH2n expected, got C{num_carbons}H{num_hydrogens}"
    
    return True, "Contains an acyclic structure with one carbon-carbon double bond fitting the formula CnH2n"