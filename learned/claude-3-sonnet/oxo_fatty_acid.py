"""
Classifies: CHEBI:59644 oxo fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is a fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"

    # Check for additional carbonyl group on the carbon chain
    chain_carbonyl_pattern = Chem.MolFromSmarts("[CX4]~[CX3](=O)~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_carbonyl_pattern):
        return False, "No additional carbonyl group on the carbon chain"

    # Check for long, unbranched carbon chain (at least 6 carbons)
    chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Carbon chain too short or branched"

    # Exclude non-fatty acid structures
    exclude_patterns = [
        Chem.MolFromSmarts("[NX3]"),  # Exclude amines
        Chem.MolFromSmarts("[OX2H]"),  # Exclude alcohols (excluding carboxylic acid)
        Chem.MolFromSmarts("[SX2]"),  # Exclude sulfur-containing groups
        Chem.MolFromSmarts("[#6;r]"),  # Exclude aromatic rings
        Chem.MolFromSmarts("[OX2][CX4][OX2]")  # Exclude sugars
    ]
    if any(mol.HasSubstructMatch(pattern) for pattern in exclude_patterns):
        return False, "Contains excluded functional groups"

    # Check oxygen count (at least 3 oxygens for carboxylic acid + carbonyl)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Too few oxygen atoms (requires at least 3)"

    return True, "Molecule is an oxo fatty acid"