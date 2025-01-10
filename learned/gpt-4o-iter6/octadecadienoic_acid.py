"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is defined as a straight-chain, C18 polyunsaturated fatty acid with two C=C double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octadecadienoic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check total carbon count allowing exact for C18 fatty acids
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 18:
        return False, f"Expected 18 carbons, found {carbon_count}"

    # Check for exactly two C=C double bonds
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE and
                            bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6)
    if double_bond_count != 2:
        return False, f"Expected 2 double bonds (C=C), found {double_bond_count}"

    # Pattern for terminal carboxylic acid group (COOH or COO-)
    terminal_carboxylic_acid_patterns = [
        Chem.MolFromSmarts("C(=O)O"),    # -COOH
        Chem.MolFromSmarts("C(=O)[O-]"), # -COO-
        Chem.MolFromSmarts("[CX3](=O)[OX1H0-1]")  # Additional for broader terminal detection
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in terminal_carboxylic_acid_patterns):
        return False, "No valid terminal carboxylic acid group found"

    return True, "Molecule is a C18 polyunsaturated fatty acid with two C=C double bonds and a terminal carboxylic acid group"