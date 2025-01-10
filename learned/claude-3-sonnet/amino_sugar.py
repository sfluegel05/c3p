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

    # More generalized sugar ring patterns
    sugar_patterns = [
        # 6-membered sugar ring (pyranose) with flexible substitution
        "[#6]1-[#6]-[#6](-[#8,#7])-[#6](-[#8,#7])-[#6](-[#8,#7])-[#6]1-[#8,#7]",
        # 5-membered sugar ring (furanose) with flexible substitution
        "[#6]1-[#6](-[#8,#7])-[#6](-[#8,#7])-[#6](-[#8,#7])-[#6]1-[#8,#7]",
        # Alternative pyranose pattern
        "[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]1",
        # Alternative furanose pattern
        "[#6]1-[#6]-[#6]-[#6]-[#6]1"
    ]

    has_sugar_ring = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sugar_ring = True
            break

    if not has_sugar_ring:
        return False, "No sugar ring structure found"

    # Enhanced amino group patterns including N-acetyl and other common modifications
    amino_patterns = [
        # Primary amine
        "[NX3;H2]",
        # Secondary amine
        "[NX3;H1][#6]",
        # N-acetyl group (very common in amino sugars)
        "[NX3;H1]C(=O)[CH3]",
        # Other N-substituted groups
        "[NX3;H1,H0]",
        # Amide nitrogen
        "[NX3][CX3](=[OX1])[#6]"
    ]

    amino_matches = []
    for pattern in amino_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat and mol.HasSubstructMatch(pat):
            amino_matches.extend(mol.GetSubstructMatches(pat))

    if not amino_matches:
        return False, "No amino groups found"

    # Check for hydroxyl groups characteristic of sugars
    hydroxyl_patterns = [
        "[OX2H1]",  # Free hydroxyl
        "[OX2]([#6])[#6]",  # Ether/acetal oxygen
    ]

    hydroxy_matches = []
    for pattern in hydroxyl_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            hydroxy_matches.extend(mol.GetSubstructMatches(pat))

    # Ring analysis
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    sugar_ring_atoms = set()
    for ring in rings:
        if len(ring) in [5, 6]:  # Furanose or pyranose
            sugar_ring_atoms.update(ring)

    # Check if amino groups are connected to the sugar ring
    amino_connected_to_ring = False
    for match in amino_matches:
        n_atom = match[0]
        for neighbor in mol.GetAtomWithIdx(n_atom).GetNeighbors():
            if neighbor.GetIdx() in sugar_ring_atoms:
                amino_connected_to_ring = True
                break
        if amino_connected_to_ring:
            break

    if not amino_connected_to_ring:
        return False, "Amino group not connected to sugar ring"

    # Final verification of composition
    ring_size = max(len(ring) for ring in rings) if rings else 0
    if ring_size not in [5, 6]:
        return False, "No suitable sugar ring size found"

    return True, "Contains amino sugar structure with amino group(s) connected to sugar ring"