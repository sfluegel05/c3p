"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: CHEBI:36195 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is any fatty acid containing at least one C=C or C#C bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for terminal carboxylic acid group
    terminal_acid_pattern = Chem.MolFromSmarts("[CH3,CH2,CH]C(=O)O")
    if not mol.HasSubstructMatch(terminal_acid_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Look for at least one C=C or C#C bond
    unsaturated_bond_pattern = Chem.MolFromSmarts("C=C,C#C")
    if not mol.HasSubstructMatch(unsaturated_bond_pattern):
        return False, "No unsaturated bonds found"

    # Check for undesirable functional groups
    exclude_patterns = [
        Chem.MolFromSmarts("O=C-O-C"),  # esters
        Chem.MolFromSmarts("C(=O)N"),   # amides
        Chem.MolFromSmarts("C1CCCCCC1") # cyclic structures
    ]
    for pattern in exclude_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains undesirable functional group"

    # Count carbons, oxygens, and check rotatable bonds
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)

    if o_count < 2 or o_count > 4:
        return False, "Incorrect number of oxygens for fatty acid"
    if n_rotatable < c_count - 2:
        return False, "Carbon chain too short for fatty acid"

    return True, "Contains at least one unsaturated bond and a terminal carboxylic acid group"