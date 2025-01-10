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
    A glycosaminoglycan is any polysaccharide containing a substantial proportion 
    of aminomonosaccharide residues.
    
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
    
    if o_count < 3:
        return False, "Too few oxygen atoms for aminomonosaccharide structure"
    
    if n_count == 0:
        return False, "No nitrogen atoms found - required for aminosugar"

    # More flexible sugar ring patterns
    sugar_patterns = [
        "[#6R1]1[#6R1][#6R1][#6R1][#6R1][#8R1]1",  # pyranose
        "[#6R1]1[#6R1][#6R1][#6R1][#8R1]1",        # furanose
        "[#6]1[#6][#6][#6][#6][#8]1",              # modified pyranose
        "[#6]1[#6][#6][#6][#8]1",                  # modified furanose
        "[#6R1]1[#6R1][#6R1][#8R1][#6R1]1"        # other sugar-like rings
    ]
    
    total_sugar_matches = 0
    for pattern in sugar_patterns:
        sugar_pattern = Chem.MolFromSmarts(pattern)
        if sugar_pattern:
            matches = len(mol.GetSubstructMatches(sugar_pattern))
            total_sugar_matches += matches

    # More flexible glycosidic linkage patterns
    glycosidic_patterns = [
        "[#6R][#8][#6R]",           # classic glycosidic
        "[#6R][#8][#6]",            # modified glycosidic
        "[#6][#8][#6R]"             # variant glycosidic
    ]
    
    glycosidic_count = 0
    for pattern in glycosidic_patterns:
        glyco_pattern = Chem.MolFromSmarts(pattern)
        if glyco_pattern:
            matches = len(mol.GetSubstructMatches(glyco_pattern))
            glycosidic_count += matches

    # Broader amino sugar patterns
    amino_sugar_patterns = [
        "[NX3,NX4][#6][#6R]1[#8][#6R][#6R][#6R][#6R]1",  # amino on/near sugar
        "[NX3,NX4][#6R]1[#8][#6R][#6R][#6R][#6R]1",      # direct amino sugar
        "[NX3,NX4]C(=O)[#6][#6R]1[#8][#6R][#6R][#6R]1",  # N-modified sugar
        "[NX3,NX4][#6][#6]1[#8][#6][#6][#6][#6]1",       # flexible amino sugar
        "[NX3,NX4][#6R]1[#8][#6][#6][#6]1"               # smaller amino sugar
    ]
    
    amino_sugar_matches = 0
    for pattern in amino_sugar_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol:
            matches = len(mol.GetSubstructMatches(pattern_mol))
            amino_sugar_matches += matches

    # Calculate ring fraction
    ring_atoms = len(Chem.GetSymmSSSR(mol))
    ring_fraction = ring_atoms / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0

    # Build features list
    features = []
    if total_sugar_matches > 0:
        features.append(f"Contains {total_sugar_matches} sugar-like rings")
    if glycosidic_count > 0:
        features.append(f"{glycosidic_count} potential glycosidic linkages")
    if amino_sugar_matches > 0:
        features.append(f"{amino_sugar_matches} amino sugar-like residues")

    # Relaxed classification criteria
    is_gag = (
        (total_sugar_matches >= 1 or glycosidic_count >= 1) and  # Sugar-like structure
        amino_sugar_matches >= 1 and                              # Must have amino sugar component
        ring_fraction >= 0.15 and                                 # Some ring content
        o_count >= 3 and                                         # Minimum oxygen
        n_count >= 1 and                                         # Must contain nitrogen
        c_count >= 6                                             # Minimum carbon
    )

    if not features:
        return False, "No characteristic glycosaminoglycan features found"
    
    reason = ("Classified as glycosaminoglycan: " + ", ".join(features)) if is_gag else \
             ("Not classified as glycosaminoglycan: " + ", ".join(features))

    return is_gag, reason