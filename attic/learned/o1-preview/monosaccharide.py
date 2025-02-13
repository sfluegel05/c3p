"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    A monosaccharide is a polyhydroxy aldehyde or ketone with three or more carbon atoms,
    and it is a single sugar unit without glycosidic connections to other units.
    This includes both open-chain and cyclic forms, as well as deoxy and amino sugars.
    
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
    
    # Check for aldehyde or ketone group in open-chain form
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    
    # Check for cyclic hemiacetal or hemiketal (cyclic form)
    hemiacetal_pattern = Chem.MolFromSmarts("[C;R][O][C;R]")
    hemiketal_pattern = Chem.MolFromSmarts("[C;R][O][C;R]")
    has_hemiacetal = mol.HasSubstructMatch(hemiacetal_pattern)
    has_hemiketal = mol.HasSubstructMatch(hemiketal_pattern)
    is_cyclic = mol.GetRingInfo().NumRings() > 0 and (has_hemiacetal or has_hemiketal)
    
    if not (has_aldehyde or has_ketone or is_cyclic):
        return False, "Molecule lacks aldehyde, ketone, or cyclic hemiacetal/hemiketal group"
    
    # Count the number of hydroxyl and amino groups attached to carbons
    hydroxyl_amino_pattern = Chem.MolFromSmarts("[C][O,N;H1,H2]")
    hydroxyl_amino_matches = mol.GetSubstructMatches(hydroxyl_amino_pattern)
    min_groups = 2  # Require at least two hydroxyl or amino groups
    if len(hydroxyl_amino_matches) < min_groups:
        return False, f"Found {len(hydroxyl_amino_matches)} hydroxyl/amino groups attached to carbons, need at least {min_groups}"
    
    # Check for glycosidic bonds to other sugar units
    # Look for anomeric carbon with exocyclic oxygen connected to another sugar
    glycosidic_bond = False
    for atom in mol.GetAtoms():
        # Identify potential anomeric carbons
        if atom.GetAtomicNum() == 6 and atom.IsInRing():
            oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if len(oxygen_neighbors) == 2:
                # One oxygen is ring oxygen, the other is exocyclic
                ring_oxygens = [
                    o for o in oxygen_neighbors
                    if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).IsInRing()
                ]
                exocyclic_oxygens = [
                    o for o in oxygen_neighbors
                    if not mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).IsInRing()
                ]
                if len(ring_oxygens) == 1 and len(exocyclic_oxygens) == 1:
                    exo_oxygen = exocyclic_oxygens[0]
                    # Check if exocyclic oxygen is connected to a non-hydrogen atom (possible glycosidic bond)
                    exo_bonds = exo_oxygen.GetBonds()
                    for bond in exo_bonds:
                        neighbor_atom = bond.GetOtherAtom(exo_oxygen)
                        if neighbor_atom.GetIdx() != atom.GetIdx() and neighbor_atom.GetAtomicNum() > 1:
                            glycosidic_bond = True
                            break
                if glycosidic_bond:
                    break
    if glycosidic_bond:
        return False, "Molecule has glycosidic linkage to another sugar unit"
    
    # If all checks passed, classify as monosaccharide
    return True, "Molecule is a monosaccharide based on structure"