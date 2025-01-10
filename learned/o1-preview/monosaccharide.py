"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or ketone with three or more carbon atoms,
    and it is a single sugar unit without glycosidic connections to other units.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove counter ions and keep only the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())
    
    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Molecule has {c_count} carbon atoms, fewer than 3"

    # Check for aldehyde or ketone group (open-chain form)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # Check for cyclic hemiacetal or hemiketal (cyclic form)
    hemiacetal_pattern = Chem.MolFromSmarts("[C;!R]=O")
    is_cyclic = mol.GetRingInfo().NumRings() > 0

    if not (has_aldehyde or has_ketone or is_cyclic):
        return False, "Molecule lacks aldehyde, ketone, or cyclic hemiacetal/hemiketal group"

    # Count the number of hydroxyl groups attached to carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[C][OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    min_hydroxyls = c_count - 2  # Monosaccharides have at least (carbon atoms - 2) hydroxyl groups
    if len(hydroxyl_matches) < min_hydroxyls:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups attached to carbons, need at least {min_hydroxyls}"
    
    # Check for glycosidic bonds to other sugar units
    # Glycosidic bonds occur when an anomeric carbon is linked to another sugar unit
    # We can check for O-glycosidic linkages that are not part of the ring
    glycosidic_pattern = Chem.MolFromSmarts("[C@H1,O]=[O,C]")  # Simplified pattern for demonstration
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings > 1:
        return False, f"Molecule has {num_rings} rings, indicating possible glycosidic linkages"
    
    # Check for multiple monosaccharide units (e.g., disaccharides)
    sugar_unit = Chem.MolFromSmarts("[C;R][O;R][C;R]")
    sugar_matches = mol.GetSubstructMatches(sugar_unit)
    if len(sugar_matches) > c_count:
        return False, "Molecule contains multiple sugar units"

    # Ensure there are no glycosidic bonds extending from the ring oxygen
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.IsInRing():
            continue
        if atom1.IsInRing() and atom2.IsInRing():
            continue
        if (atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8):  # Oxygen atom
            # Check if the oxygen is connected to carbons outside the ring
            if (atom1.IsInRing() and not atom2.IsInRing()) or (atom2.IsInRing() and not atom1.IsInRing()):
                return False, "Molecule has glycosidic linkage to another unit"
    
    return True, "Molecule is a monosaccharide based on structure"