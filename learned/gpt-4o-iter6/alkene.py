"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is defined as an acyclic branched or unbranched hydrocarbon having one C=C double bond,
    following the general formula CnH2n.

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
    
    # Check for cycles; alkenes must be acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; alkenes must be acyclic"
    
    # Count carbon-carbon double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() 
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and 
                            bond.GetBeginAtom().GetSymbol() == 'C' and 
                            bond.GetEndAtom().GetSymbol() == 'C')
    
    # Alkenes must have exactly one C=C double bond
    if double_bond_count != 1:
        return False, f"Expected 1 C=C double bond, found {double_bond_count}"
    
    # Verify the formula CnH2n
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    
    num_carbons = atom_counts.get('C', 0)
    num_hydrogens = atom_counts.get('H', 0)
    
    if num_hydrogens != 2 * num_carbons:
        return False, f"Formula check failed: CnH2n expected, got C{num_carbons}H{num_hydrogens}"
    
    return True, "Contains an acyclic structure with one carbon-carbon double bond fitting the formula CnH2n"