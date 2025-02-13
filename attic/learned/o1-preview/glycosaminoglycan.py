"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is any polysaccharide containing a substantial proportion of aminomonosaccharide residues.
    
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

    # Define SMARTS patterns for 5- and 6-membered sugar rings
    sugar_patterns = [
        Chem.MolFromSmarts("[C;R]1([O;R][C;R][C;R][C;R][C;R]1)"),  # 6-membered ring
        Chem.MolFromSmarts("[C;R]1([O;R][C;R][C;R][C;R]1)")        # 5-membered ring
    ]
    # Define SMARTS pattern for amino group attached to ring carbon
    amino_group = Chem.MolFromSmarts("[C;R][N;!R]")

    # Check for valid SMARTS patterns
    if None in sugar_patterns or amino_group is None:
        return False, "Error in SMARTS patterns"

    # Find sugar rings
    sugar_ring_matches = []
    for pattern in sugar_patterns:
        matches = mol.GetSubstructMatches(pattern)
        sugar_ring_matches.extend(matches)

    if len(sugar_ring_matches) == 0:
        return False, "No sugar rings found"

    # Check for amino sugars
    contains_amino_sugar = False
    for match in sugar_ring_matches:
        ring_atoms = set(match)
        # Search for amino group attached to ring
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7 and neighbor.GetIdx() not in ring_atoms:
                    contains_amino_sugar = True
                    break
            if contains_amino_sugar:
                break
        if contains_amino_sugar:
            break

    if not contains_amino_sugar:
        return False, "No amino sugars found"

    # Check for polysaccharide (multiple sugar rings connected via glycosidic bonds)
    if len(sugar_ring_matches) < 2:
        return False, "Not enough sugar rings for polysaccharide"

    # Check for glycosidic bonds between sugar rings (oxygen bridges)
    glycosidic_bonds = 0
    for i in range(len(sugar_ring_matches)):
        for j in range(i+1, len(sugar_ring_matches)):
            ring1 = set(sugar_ring_matches[i])
            ring2 = set(sugar_ring_matches[j])
            # Find connecting atoms
            for atom_idx1 in ring1:
                atom1 = mol.GetAtomWithIdx(atom_idx1)
                for neighbor in atom1.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in ring1:
                        # Neighbor is oxygen not in ring1
                        for neighbor2 in neighbor.GetNeighbors():
                            if neighbor2.GetIdx() in ring2:
                                glycosidic_bonds += 1
                                break
                if glycosidic_bonds > 0:
                    break
            if glycosidic_bonds > 0:
                break
        if glycosidic_bonds > 0:
            break

    if glycosidic_bonds == 0:
        return False, "No glycosidic bonds between sugar rings"

    return True, "Contains polysaccharide with amino sugar residues"