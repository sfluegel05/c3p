"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: CHEBI:32395 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid has one double or triple bond in the carbon chain and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count double and triple bonds in the carbon chain, excluding the carboxylic acid group
    double_bonds = 0
    triple_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Exclude the double bond in the carboxylic acid group
            if not (bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 8):
                double_bonds += 1
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            triple_bonds += 1

    total_unsaturations = double_bonds + triple_bonds

    # Check for exactly one double or triple bond in the carbon chain
    if total_unsaturations != 1:
        return False, f"Found {total_unsaturations} unsaturations in the carbon chain, need exactly 1"

    # Ensure the unsaturation is in the carbon chain, not in a ring
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            if bond.GetBeginAtom().IsInRing() or bond.GetEndAtom().IsInRing():
                return False, "Unsaturation is in a ring, not in the carbon chain"

    # Check that the molecule is a simple chain (no rings, no branching beyond the carboxylic acid)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a simple fatty acid chain"

    # Check that the molecule is a simple chain (no branching beyond the carboxylic acid)
    # Count the number of carbons with more than 2 connections
    branching_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2:
            branching_carbons += 1
    if branching_carbons > 1:  # Allow one branching point for the carboxylic acid
        return False, "Molecule is too branched to be a simple fatty acid"

    # Check for additional functional groups that would disqualify the molecule
    # Allow hydroxyl groups, but no other functional groups
    allowed_functional_groups = ["[OH]", "[OX2H1]"]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1:  # Oxygen with one bond (hydroxyl)
            continue
        elif atom.GetAtomicNum() == 6 and atom.GetDegree() <= 2:  # Carbon in the chain
            continue
        else:
            return False, "Molecule contains additional functional groups"

    return True, "Contains a carboxylic acid group and exactly one double or triple bond in the carbon chain"