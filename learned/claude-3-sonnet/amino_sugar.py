"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:22044 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar having one or more alcoholic hydroxy groups 
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible sugar ring patterns that allow for substitutions
    pyranose_pattern = Chem.MolFromSmarts("[CR0]1[CR0][CR0]([OR0,NR0;!H0,$(NC=O)])[CR0]([OR0,NR0;!H0,$(NC=O)])[CR0]([OR0,NR0;!H0,$(NC=O)])[CR0]1[OR0,NR0;!H0,$(NC=O)]")
    furanose_pattern = Chem.MolFromSmarts("[CR0]1[CR0]([OR0,NR0;!H0,$(NC=O)])[CR0]([OR0,NR0;!H0,$(NC=O)])[CR0]([OR0,NR0;!H0,$(NC=O)])[CR0]1[OR0,NR0;!H0,$(NC=O)]")
    
    # Alternative sugar patterns for complex cases
    alt_sugar_pattern = Chem.MolFromSmarts("[CR0]1[CR0][CR0]([OR0,NR0])[CR0]([OR0,NR0])[CR0]([OR0,NR0])[OR0,NR0]1")
    
    has_sugar_ring = any([
        mol.HasSubstructMatch(pyranose_pattern),
        mol.HasSubstructMatch(furanose_pattern),
        mol.HasSubstructMatch(alt_sugar_pattern)
    ])
    
    if not has_sugar_ring:
        return False, "No sugar ring structure found"

    # Enhanced amino group patterns
    amino_patterns = [
        Chem.MolFromSmarts("[NX3;H2,H1]"),  # Primary/secondary amines
        Chem.MolFromSmarts("[NX3;H1;$(NC=O)]"),  # Amides including N-acetyl
        Chem.MolFromSmarts("[NX3;H0;$(N(C)C=O)]"),  # N-substituted amides
        Chem.MolFromSmarts("[NX3;$(NC([#6])[#6])]")  # Other N-substituted groups
    ]
    
    has_amino = any(mol.HasSubstructMatch(pattern) for pattern in amino_patterns)
    
    if not has_amino:
        return False, "No amino or substituted amino groups found"

    # Check for characteristic sugar features
    sugar_features = [
        ("[CR0]([OR0,NR0])([OR0,NR0])", "sugar carbon with multiple heteroatom substituents"),
        ("[CR0]1[OR0][CR0][CR0][CR0]([OR0,NR0])[CR0]1", "pyranose ring with heteroatom substituent"),
        ("[CR0]([OR0])([CR0])[CR0]([OR0,NR0])", "sugar carbon chain with heteroatom substituents")
    ]
    
    found_features = []
    for pattern, feature in sugar_features:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_features.append(feature)
    
    if not found_features:
        return False, "Missing characteristic sugar structural features"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 1:
        return False, "Too few hydroxyl groups for a sugar"

    # Verify basic composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 3 or o_count < 1 or n_count < 1:
        return False, "Insufficient atoms for amino sugar structure"

    # Additional check for cyclic structure
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"

    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size in [5,6] for size in ring_sizes):
        return False, "No suitable sugar ring size found"

    return True, f"Contains amino sugar features: {', '.join(found_features)}"