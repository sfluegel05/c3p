"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is defined as an acyclic branched or unbranched hydrocarbon having 
    exactly one carbon-carbon double bond with the general formula CnH2n.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains cycles"

    # Look for carbon-carbon double bonds (C=C)
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    num_double_bonds = len(double_bond_matches)

    if num_double_bonds != 1:
        return False, f"Found {num_double_bonds} carbon-carbon double bonds, need exactly 1"

    # Check if the molecular formula matches CnH2n
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)

    if h_count != 2 * c_count:
        return False, f"Molecular formula does not match CnH2n (found C{c_count}H{h_count})"

    return True, "Molecule is an alkene with one carbon-carbon double bond and fits the general formula CnH2n"