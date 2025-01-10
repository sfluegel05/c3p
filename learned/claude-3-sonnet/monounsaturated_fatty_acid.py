"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    Must have exactly one double/triple bond and a carboxylic acid group.

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

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Count number of carboxylic acid groups
    carboxylic_matches = len(mol.GetSubstructMatches(carboxylic_pattern))
    if carboxylic_matches > 1:
        return False, "Multiple carboxylic acid groups found"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = len(mol.GetSubstructMatches(double_bond_pattern))
    
    # Count triple bonds
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    triple_bonds = len(mol.GetSubstructMatches(triple_bond_pattern))
    
    # Total unsaturations (excluding the carboxylic acid)
    total_unsaturations = double_bonds + triple_bonds
    
    if total_unsaturations == 0:
        return False, "No double or triple bonds found in carbon chain"
    elif total_unsaturations > 1:
        return False, f"Found {total_unsaturations} unsaturations, expected 1"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty acid"
    
    # Check for cyclic structures (excluding small rings that might appear in some MUFAs)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if any(size > 4 for size in ring_sizes):
        return False, "Contains large ring structure"

    # Count branching points
    branching_pattern = Chem.MolFromSmarts("[*]([*])([*])([*])[*]")
    branching_points = len(mol.GetSubstructMatches(branching_pattern))
    if branching_points > 1:  # Allow some branching for methyl groups
        return False, "Too many branching points for fatty acid"

    # Additional check for conjugated systems
    conjugated_pattern = Chem.MolFromSmarts("C=C-C=C")
    if mol.HasSubstructMatch(conjugated_pattern):
        return False, "Contains conjugated double bonds"

    return True, "Contains one unsaturation and carboxylic acid group in appropriate structure"