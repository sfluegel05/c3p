"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: CHEBI:17855 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides are N-acyl-sphingoid bases with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingoid base backbone patterns
    sphingoid_patterns = [
        Chem.MolFromSmarts("[N;X3][C;X4][C;X4][C;X4][C;X4][C;X4]"),  # Linear
        Chem.MolFromSmarts("[N;X3][C;X4][C;X4][C;X4][C;X4]=[C;X3]"),  # With double bond
        Chem.MolFromSmarts("[N;X3][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]")  # Longer chain
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in sphingoid_patterns):
        return False, "No sphingoid base backbone found"

    # Look for amide group
    amide_pattern = Chem.MolFromSmarts("[N;X3][C;X3](=[O;X1])")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Look for fatty acid chains
    fatty_acid_patterns = [
        Chem.MolFromSmarts("[C;X4][C;X4][C;X4][C;X4]"),  # Minimum length
        Chem.MolFromSmarts("[C;X4][C;X4][C;X4][C;X4][C;X4]"),  # Longer chain
        Chem.MolFromSmarts("[C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]"),  # Even longer chain
        Chem.MolFromSmarts("[C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]")  # Longest chain
    ]
    fatty_acid_matches = sum(mol.GetSubstructMatches(pat) for pat in fatty_acid_patterns)
    if fatty_acid_matches < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds as a supplementary check
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"

    # Count carbons, nitrogens, and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 16:
        return False, "Too few carbons for ceramide"
    if n_count != 2:
        return False, "Must have exactly 2 nitrogens (1 amide, 1 sphingoid base)"
    if o_count < 3:
        return False, "Too few oxygens for ceramide"

    return True, "Contains sphingoid base with fatty acid chain attached via amide bond"