"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: glucosinolate compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N-O-SO3 group (both charged and uncharged forms)
    # More flexible pattern to catch variations
    sulfate_patterns = [
        Chem.MolFromSmarts('[NX2]-[OX2]-[SX4](=[OX1])(=[OX1])[OX1,OX2H1,OX2-]'),
        Chem.MolFromSmarts('[NX2]-[OX2]-[SX4]([OX1,OX2H1,OX2-])(=[OX1])=[OX1]')
    ]
    has_sulfate = any(mol.HasSubstructMatch(pattern) for pattern in sulfate_patterns)
    if not has_sulfate:
        return False, "Missing N-O-SO3 group"

    # Look for C=N-O linkage with some flexibility
    cn_bond_pattern = Chem.MolFromSmarts('[CX3]=[NX2]-[OX2]')
    if not mol.HasSubstructMatch(cn_bond_pattern):
        return False, "Missing C=N-O linkage"

    # Check for thioglycosidic linkage (C-S-C=N)
    # More flexible pattern to catch variations
    thio_patterns = [
        Chem.MolFromSmarts('[CX4]-[SX2]-[CX3]=[NX2]'),
        Chem.MolFromSmarts('[#6]-[SX2]-[CX3]=[NX2]')
    ]
    has_thio = any(mol.HasSubstructMatch(pattern) for pattern in thio_patterns)
    if not has_thio:
        return False, "Missing thioglycosidic linkage"

    # Check for glucose moiety with more flexible patterns
    sugar_patterns = [
        # Pattern for beta-D-glucopyranose with various representations
        Chem.MolFromSmarts('[OX2][#6]1[OX2][#6]([SX2])[#6]([OX2])[#6]([OX2])[#6]([OX2])[#6]1[CX4][OX2]'),
        Chem.MolFromSmarts('[OX2][#6]-1-[OX2]-[#6](-[SX2])-[#6](-[OX2])-[#6](-[OX2])-[#6](-[OX2])-[#6]-1-[CX4]-[OX2]')
    ]
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "Missing glucose moiety"

    # Count key atoms to ensure reasonable composition
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if s_count < 2:  # Need at least 2 sulfur atoms (thioglycosidic + sulfate)
        return False, "Insufficient sulfur atoms"
    if o_count < 7:  # Need multiple oxygens for sugar + sulfate
        return False, "Insufficient oxygen atoms"
    if n_count != 1:  # Need exactly 1 nitrogen
        return False, "Incorrect number of nitrogen atoms"

    # Verify the overall connectivity but with more flexible pattern
    core_pattern = Chem.MolFromSmarts('[#6]-[SX2]-[CX3](=[NX2]-[OX2]-[SX4]~[OX1,OX2H1,OX2-])-[#6]')
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Incorrect overall connectivity"

    return True, "Contains complete glucosinolate structure with correct connectivity"