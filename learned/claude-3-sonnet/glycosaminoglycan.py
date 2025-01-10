"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count key atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if o_count < 3:  # Relaxed oxygen requirement
        return False, "Too few oxygen atoms for aminomonosaccharide structure"
    
    if n_count == 0:  # Still need nitrogen for amino groups
        return False, "No nitrogen atoms found - required for aminosugar"

    # Look for various sugar-like ring patterns
    patterns = [
        "[CR1]1[CR1][CR1][CR1][CR1]O1",  # pyranose
        "[CR1]1[CR1][CR1][CR1]O1",       # furanose
        "[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]O1",  # modified pyranose
        "[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1]O1"            # modified furanose
    ]
    
    total_sugar_matches = 0
    for pattern in patterns:
        sugar_pattern = Chem.MolFromSmarts(pattern)
        if sugar_pattern:
            matches = len(mol.GetSubstructMatches(sugar_pattern))
            total_sugar_matches += matches

    # Look for amino groups in various contexts
    amino_patterns = [
        "[NX3H,NX3H2][CH1,CH2][OH1,OR]",  # classic amino sugar
        "[NX3H,NX3H2][CR0,CR1][OR0,OR1]",  # modified amino sugar
        "[NX3H,NX3H2]C(=O)",               # N-acetyl group
        "[NX3]C[CR1]1O[CR1][CR1][CR1][CR1]1"  # N-substituted sugar
    ]
    
    total_amino_matches = 0
    for pattern in amino_patterns:
        amino_pattern = Chem.MolFromSmarts(pattern)
        if amino_pattern:
            matches = len(mol.GetSubstructMatches(amino_pattern))
            total_amino_matches += matches

    # Look for characteristic modifications
    sulfate = Chem.MolFromSmarts("OS(=O)(=O)[OH1,O-]")
    carboxyl = Chem.MolFromSmarts("C(=O)[OH1,O-]")
    acetyl = Chem.MolFromSmarts("NC(=O)C")
    
    has_sulfate = mol.HasSubstructMatch(sulfate) if sulfate else False
    has_carboxyl = mol.HasSubstructMatch(carboxyl) if carboxyl else False
    has_acetyl = mol.HasSubstructMatch(acetyl) if acetyl else False

    # Build features list
    features = []
    if total_sugar_matches > 0:
        features.append(f"Contains {total_sugar_matches} sugar-like rings")
    if total_amino_matches > 0:
        features.append(f"{total_amino_matches} amino groups")
    if has_sulfate:
        features.append("sulfate groups")
    if has_carboxyl:
        features.append("carboxyl groups")
    if has_acetyl:
        features.append("acetyl groups")

    # Classification criteria
    is_gag = (
        (total_sugar_matches > 0 or (o_count >= 4 and n_count >= 1)) and  # Either has sugar rings or enough O/N atoms
        total_amino_matches > 0 and                                        # Must have amino groups
        (has_sulfate or has_carboxyl or has_acetyl)                       # Must have at least one characteristic modification
    )

    if not features:
        return False, "No characteristic glycosaminoglycan features found"
    
    reason = ("Classified as glycosaminoglycan: " + ", ".join(features)) if is_gag else \
             ("Not classified as glycosaminoglycan: " + ", ".join(features))

    return is_gag, reason