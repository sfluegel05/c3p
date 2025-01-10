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

    # Look for pyranose (6-membered) or furanose (5-membered) sugar ring patterns
    # More flexible pattern that matches sugar-like rings with various substituents
    pyranose_pattern = Chem.MolFromSmarts("[C]1[C][C]([O,N;!H0,$(NC=O)])[C]([O,N;!H0,$(NC=O)])[C]([O,N;!H0,$(NC=O)])[C]1[O,N;!H0,$(NC=O)]")
    furanose_pattern = Chem.MolFromSmarts("[C]1[C]([O,N;!H0,$(NC=O)])[C]([O,N;!H0,$(NC=O)])[C]([O,N;!H0,$(NC=O)])[C]1[O,N;!H0,$(NC=O)]")
    
    has_sugar_ring = mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)
    if not has_sugar_ring:
        return False, "No sugar ring structure found"

    # Look for amino groups, including N-acetyl groups which are common in amino sugars
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Primary/secondary amines
    n_acetyl_pattern = Chem.MolFromSmarts("[NX3;H1;$(NC(=O)C)]")  # N-acetyl groups
    amide_pattern = Chem.MolFromSmarts("[NX3;H1;$(NC=O)]")  # Other amides
    
    has_amine = mol.HasSubstructMatch(amine_pattern)
    has_n_acetyl = mol.HasSubstructMatch(n_acetyl_pattern)
    has_amide = mol.HasSubstructMatch(amide_pattern)
    
    if not (has_amine or has_n_acetyl or has_amide):
        return False, "No amino or substituted amino groups found"

    # Count hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OH1]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    if hydroxyl_matches < 2:
        return False, "Too few hydroxyl groups for a sugar"

    # Verify cyclic nature
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"
    
    # Check ring sizes - sugars typically have 5 or 6-membered rings
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size in [5,6] for size in ring_sizes):
        return False, "No suitable sugar ring size found"

    # Count carbons and oxygens in ring systems
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.IsInRing())
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.IsInRing())
    
    if c_count < 5:
        return False, "Too few carbons in ring system"
    if o_count < 1:
        return False, "Too few oxygens in ring system"

    # Check for characteristic sugar carbon with attached oxygen/nitrogen
    sugar_carbon_pattern = Chem.MolFromSmarts("[C;R]([O,N])")
    if not mol.HasSubstructMatch(sugar_carbon_pattern):
        return False, "Missing characteristic sugar carbon-heteroatom bond"

    return True, "Contains sugar ring structure with amino group(s) replacing hydroxyl position(s)"