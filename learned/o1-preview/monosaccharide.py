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
    This includes both open-chain and cyclic forms.
    
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
    is_cyclic = mol.GetRingInfo().NumRings() > 0
    
    if not (has_aldehyde or has_ketone or is_cyclic):
        return False, "Molecule lacks aldehyde, ketone, or cyclic hemiacetal/hemiketal group"
    
    # Count the number of hydroxyl groups attached to carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[C][O;H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    min_hydroxyls = c_count - 2  # Monosaccharides have at least (carbon atoms - 2) hydroxyl groups
    if len(hydroxyl_matches) < min_hydroxyls:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups attached to carbons, need at least {min_hydroxyls}"
    
    # Check for glycosidic bonds to other sugar units
    # Analyze anomeric carbons to detect glycosidic linkages
    glycosidic_bond = False
    for atom in mol.GetAtoms():
        # Identify potential anomeric carbons
        if atom.GetAtomicNum() == 6 and atom.IsInRing():
            oxygen_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if len(oxygen_neighbors) == 2:
                # Check if one oxygen is in ring (ring oxygen) and the other is exocyclic
                ring_oxygens = [o for o in oxygen_neighbors if mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).IsInRing()]
                exocyclic_oxygens = [o for o in oxygen_neighbors if not mol.GetBondBetweenAtoms(atom.GetIdx(), o.GetIdx()).IsInRing()]
                if len(ring_oxygens) == 1 and len(exocyclic_oxygens) == 1:
                    exo_oxygen = exocyclic_oxygens[0]
                    # Check if exocyclic oxygen is connected to hydrogen (OH group)
                    exo_oxygen_neighbors = [nbr for nbr in exo_oxygen.GetNeighbors() if nbr.GetIdx() != atom.GetIdx()]
                    if len(exo_oxygen_neighbors) == 0:
                        # Oxygen has no other heavy atom neighbors, likely missing explicit hydrogen
                        num_hydrogens = exo_oxygen.GetTotalNumHs()
                        if num_hydrogens == 1:
                            continue  # OH group, acceptable
                        else:
                            return False, "Exocyclic oxygen on anomeric carbon lacks hydrogen"
                    else:
                        # Exocyclic oxygen is connected to another atom
                        for nbr in exo_oxygen_neighbors:
                            if nbr.GetAtomicNum() != 1:
                                # Oxygen is connected to a non-hydrogen atom, indicating glycosidic linkage
                                glycosidic_bond = True
                                break
                if glycosidic_bond:
                    break
    if glycosidic_bond:
        return False, "Molecule has glycosidic linkage to another sugar unit"
    
    # If all checks passed, classify as monosaccharide
    return True, "Molecule is a monosaccharide based on structure"