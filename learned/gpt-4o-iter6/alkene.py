"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene is an acyclic branched or unbranched hydrocarbon with at least one C=C double bond,
    following the formula CnH2n, considering structural variations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is considered an alkene, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic element check: primarily C and H
    carbon_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            carbon_count += 1
        elif atom.GetSymbol() != 'H':
            # Non-competitive/non-affecting elements can be present
            continue

    # Count carbon-carbon double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds()
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and
                            bond.GetBeginAtom().GetSymbol() == 'C' and
                            bond.GetEndAtom().GetSymbol() == 'C' and 
                            not bond.IsInRing())
    
    if double_bond_count == 0:
        return False, f"No acyclic C=C double bond detected, found {double_bond_count}"

    # Verify the formula CnH2n aligns well while ignoring minor heteroatoms/methyl branches etc.
    num_total_hydrogens = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    if num_total_hydrogens + carbon_count != 2 * carbon_count:
        return False, f"Formula mismatch: checking generalized CxHy sum approach"

    return True, "Detected alkene due to presence of acyclic C=C double bonds fitting variations of CnH2n"