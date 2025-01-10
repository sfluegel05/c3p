"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Verify the presence of 18 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Expected 18 carbon atoms, found {c_count}"

    # Check for exactly two C=C double bonds in the structure
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE and
                            bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6)
    if double_bond_count != 2:
        return False, f"Expected 2 double bonds (C=C), found {double_bond_count}"

    # Check for terminal carboxylic acid group or its anionic form
    carboxylic_group_patterns = [
        Chem.MolFromSmarts("C(=O)[OX1H1]"),    # -COOH
        Chem.MolFromSmarts("C(=O)[O-]"),       # -COO-
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in carboxylic_group_patterns):
        return False, "No valid terminal carboxylic acid group found"

    return True, "Molecule is a C18 polyunsaturated fatty acid with two C=C double bonds and a terminal carboxylic acid group"