"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: diketone
A compound containing exactly two ketone functionalities
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone based on its SMILES string.
    A diketone contains exactly two ketone groups (C=O where C is bonded to carbons).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all carbonyls
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])")
    if carbonyl_pattern is None:
        return False, "Error in SMARTS pattern"
    
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 2:
        return False, f"Found only {len(carbonyl_matches)} carbonyl groups"

    # Find carboxylic acids
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    acid_matches = set() if acid_pattern is None else set(x[0] for x in mol.GetSubstructMatches(acid_pattern))
    
    # Find esters
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
    ester_matches = set() if ester_pattern is None else set(x[0] for x in mol.GetSubstructMatches(ester_pattern))
    
    # Find amides
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    amide_matches = set() if amide_pattern is None else set(x[0] for x in mol.GetSubstructMatches(amide_pattern))

    # Get all non-ketone carbonyls
    non_ketone_carbonyls = acid_matches | ester_matches | amide_matches
    
    # Get ketone carbonyls (carbonyls that aren't part of acids/esters/amides)
    ketone_carbons = set(match[0] for match in carbonyl_matches) - non_ketone_carbonyls
    
    if len(ketone_carbons) != 2:
        return False, f"Found {len(ketone_carbons)} ketone groups, need exactly 2"

    # Validate each ketone carbon has at least one carbon neighbor
    for carbon_idx in ketone_carbons:
        atom = mol.GetAtomWithIdx(carbon_idx)
        carbon_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
        if carbon_neighbors == 0:
            return False, "Ketone carbon must be bonded to at least one carbon atom"

    # Check for specific diketone patterns
    alpha_diketone = Chem.MolFromSmarts("[#6]C(=O)C(=O)[#6]")
    cyclic_13_diketone = Chem.MolFromSmarts("[#6]1[#6]C(=O)[#6][#6]C1=O")
    cyclic_14_diketone = Chem.MolFromSmarts("[#6]1[#6]C(=O)[#6][#6]C(=O)[#6]1")
    
    pattern_found = ""
    if mol.HasSubstructMatch(alpha_diketone):
        pattern_found = " (alpha-diketone pattern)"
    elif mol.HasSubstructMatch(cyclic_13_diketone):
        pattern_found = " (1,3-diketone pattern)"
    elif mol.HasSubstructMatch(cyclic_14_diketone):
        pattern_found = " (1,4-diketone pattern)"

    return True, f"Contains exactly two ketone groups{pattern_found}"