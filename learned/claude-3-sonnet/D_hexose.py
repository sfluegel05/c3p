"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:16234 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    D-hexose is a hexose (6-carbon monosaccharide) with D-configuration at C5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms to ensure basic composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count != 6:
        return False, f"Must have exactly 6 carbons, found {c_count}"
    if o_count != 6:
        return False, f"Must have exactly 6 oxygens, found {o_count}"
    
    # Check for any non C/H/O atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, "Contains elements other than C, H, O"

    # Check for carboxyl groups (would indicate uronic acid)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    if mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Contains carboxyl group (uronic acid)"

    # Check for non-hydroxyl oxygen (ethers, except for hemiacetal)
    ether_pattern = Chem.MolFromSmarts("[OX2]([CX4])[CX4]")
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    ring_o_pattern = Chem.MolFromSmarts("[OR0]")
    ring_o_matches = len(mol.GetSubstructMatches(ring_o_pattern))
    
    if ether_matches > ring_o_matches:
        return False, "Contains non-hemiacetal ether linkages"

    # Patterns for different forms of D-hexoses
    patterns = [
        # Open chain D-aldohexose
        "[H]C(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",  # D-glucose
        "[H]C(=O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)CO",  # D-galactose
        "[H]C(=O)[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)CO", # D-mannose
        "[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO", # D-allose
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO", # D-altrose
        "[H]C(=O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)CO",  # D-gulose
        "[H]C(=O)[C@H](O)[C@H](O)[C@H](O)[C@H](O)CO",   # D-idose
        "[H]C(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",  # D-talose
        
        # α-D-pyranose forms (hemiacetal OH is axial)
        "O[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",  # α-D-glucopyranose
        "O[C@H]1O[C@H](CO)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O",  # α-D-galactopyranose
        "O[C@H]1O[C@H](CO)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O", # α-D-mannopyranose
        
        # β-D-pyranose forms (hemiacetal OH is equatorial)
        "O[C@H]1O[C@@H](CO)[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",  # β-D-glucopyranose
        "O[C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O",  # β-D-galactopyranose
        "O[C@H]1O[C@@H](CO)[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O", # β-D-mannopyranose
        
        # D-furanose forms
        "O1[C@@H]([C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CO",  # D-glucofuranose
        "O1[C@H]([C@H](O)[C@H](O)[C@@H]1O)[C@H](O)CO",   # D-galactofuranose
    ]

    # Convert patterns to molecules
    pattern_mols = [Chem.MolFromSmarts(p) for p in patterns]
    
    # Check if molecule matches any pattern
    for pattern in pattern_mols:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, "Matches D-hexose pattern with correct stereochemistry"

    # If no pattern matched but basic composition is correct,
    # it might be a valid D-hexose in a different conformation
    chiral_centers = Chem.FindMolChiralCenters(mol)
    if len(chiral_centers) >= 4:
        return False, "Has correct composition but stereochemistry doesn't match D-hexose patterns"

    return False, "Does not match D-hexose pattern"