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
    A glycosaminoglycan is a polysaccharide containing aminomonosaccharide residues.
    
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
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if o_count < 4:
        return False, "Too few oxygen atoms for polysaccharide structure"
    
    if n_count == 0:
        return False, "No nitrogen atoms found - required for aminosugar"

    # Look for sugar ring patterns
    sugar_patterns = [
        "[CR1]1[CR1][CR1][CR1][CR1]O1",  # pyranose
        "[CR1]1[CR1][CR1][CR1]O1",       # furanose
        "[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]O1",  # modified pyranose
        "[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1]O1"            # modified furanose
    ]
    
    total_sugar_matches = 0
    for pattern in sugar_patterns:
        sugar_pattern = Chem.MolFromSmarts(pattern)
        if sugar_pattern:
            matches = len(mol.GetSubstructMatches(sugar_pattern))
            total_sugar_matches += matches

    # Look for glycosidic linkages
    glycosidic = Chem.MolFromSmarts("[CR1]O[CR1]")
    glycosidic_count = len(mol.GetSubstructMatches(glycosidic)) if glycosidic else 0

    # Look for amino sugar patterns
    amino_sugar_patterns = [
        "[NX3H,NX3H2][CH1,CH2][CR1]1O[CR1][CR1][CR1][CR1]1",  # amino directly on sugar
        "[NX3H,NX3H2][CR0,CR1]1O[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]1",  # modified amino sugar
        "[NX3]C(=O)C[CR1]1O[CR1][CR1][CR1][CR1]1"  # N-acetylated sugar
    ]
    
    amino_sugar_matches = 0
    for pattern in amino_sugar_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol:
            matches = len(mol.GetSubstructMatches(pattern_mol))
            amino_sugar_matches += matches

    # Look for characteristic modifications
    sulfate = Chem.MolFromSmarts("OS(=O)(=O)[OH1,O-]")
    carboxyl = Chem.MolFromSmarts("C(=O)[OH1,O-]")
    acetyl = Chem.MolFromSmarts("NC(=O)C")
    
    has_sulfate = mol.HasSubstructMatch(sulfate) if sulfate else False
    has_carboxyl = mol.HasSubstructMatch(carboxyl) if carboxyl else False
    has_acetyl = mol.HasSubstructMatch(acetyl) if acetyl else False

    # Calculate ring fraction to distinguish from peptides
    ring_atoms = len(Chem.GetSymmSSSR(mol))
    ring_fraction = ring_atoms / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0

    # Build features list
    features = []
    if total_sugar_matches > 0:
        features.append(f"Contains {total_sugar_matches} sugar rings")
    if glycosidic_count > 0:
        features.append(f"{glycosidic_count} glycosidic linkages")
    if amino_sugar_matches > 0:
        features.append(f"{amino_sugar_matches} amino sugar residues")
    if has_sulfate:
        features.append("sulfate groups")
    if has_carboxyl:
        features.append("carboxyl groups")
    if has_acetyl:
        features.append("acetyl groups")

    # Classification criteria - more stringent requirements
    is_gag = (
        total_sugar_matches >= 1 and                    # Must have at least one sugar ring
        glycosidic_count >= 1 and                      # Must have glycosidic linkages
        amino_sugar_matches >= 1 and                    # Must have amino sugar residues
        (has_sulfate or has_carboxyl or has_acetyl) and # Must have characteristic modifications
        ring_fraction >= 0.2 and                        # Significant ring content
        o_count >= 4 and                                # Minimum oxygen requirement
        c_count >= 12                                   # Minimum carbon requirement
    )

    if not features:
        return False, "No characteristic glycosaminoglycan features found"
    
    reason = ("Classified as glycosaminoglycan: " + ", ".join(features)) if is_gag else \
             ("Not classified as glycosaminoglycan: " + ", ".join(features))

    return is_gag, reason