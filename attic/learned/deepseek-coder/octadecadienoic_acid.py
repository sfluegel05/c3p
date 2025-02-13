"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: CHEBI:32395 octadecadienoic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_octadecadienoic_acid(smiles: str):
    """
    Determines if a molecule is an octadecadienoic acid based on its SMILES string.
    An octadecadienoic acid is a straight-chain, C18 polyunsaturated fatty acid with two C=C double bonds.

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

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 18:
        return False, f"Expected 18 carbon atoms, found {c_count}"

    # Count the number of C=C double bonds (excluding carbonyls)
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Check if the double bond is between two carbons (C=C)
            if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
                double_bond_count += 1

    if double_bond_count != 2:
        return False, f"Expected 2 C=C double bonds, found {double_bond_count}"

    # Check for disallowed functional groups (e.g., hydroperoxides, triple bonds, rings)
    disallowed_patterns = [
        Chem.MolFromSmarts("[OX2][OX2]"),  # Hydroperoxides
        Chem.MolFromSmarts("[CX2]#[CX2]"),  # Triple bonds
        Chem.MolFromSmarts("[R]"),  # Rings
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Molecule contains disallowed functional groups or rings"

    # Check if the molecule is a straight-chain fatty acid with minimal branching
    # Allow up to 2 branching carbons (e.g., for stereochemistry or slight branching)
    branching_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if atom.GetDegree() > 2:
                branching_carbons += 1

    if branching_carbons > 2:
        return False, "Molecule has too much branching to be considered a straight-chain fatty acid"

    return True, "C18 fatty acid with 2 C=C double bonds, a carboxylic acid group, and minimal branching"