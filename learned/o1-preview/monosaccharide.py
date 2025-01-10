"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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
    
    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Molecule has {c_count} carbon atoms, fewer than 3"

    # Check for glycosidic bonds (ether linkages between two carbons)
    # Glycosidic bond pattern: C-O-C where both carbons are anomeric centers
    # Since monosaccharides should not have glycosidic linkages, we can check for C-O-C bonds
    glycosidic_pattern = Chem.MolFromSmarts("[C]-O-[C]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if glycosidic_matches:
        return False, "Molecule has glycosidic linkages, not a single sugar unit"

    # Check for aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

    # Check for ketone group
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    # Check for hemiacetal or hemiketal (in cyclic forms)
    hemiacetal_pattern = Chem.MolFromSmarts("[C;H1,H2](O)[O][C]")
    hemiacetal_matches = mol.GetSubstructMatches(hemiacetal_pattern)

    # Check if molecule has either aldehyde, ketone, or hemiacetal group
    if not (aldehyde_matches or ketone_matches or hemiacetal_matches):
        return False, "No aldehyde, ketone, or hemiacetal group found"

    # Check for multiple hydroxyl groups attached to carbons
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups attached to carbons, need at least 2"

    # Check for cyclic forms (furanose and pyranose rings)
    # Furanose ring pattern: 5-membered ring with 4 carbons and 1 oxygen
    furanose_pattern = Chem.MolFromSmarts("C1COC(C1)O")
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    # Pyranose ring pattern: 6-membered ring with 5 carbons and 1 oxygen
    pyranose_pattern = Chem.MolFromSmarts("C1CCOC(C1)O")
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    if not (furanose_matches or pyranose_matches):
        # If not cyclic, check for open-chain form with multiple hydroxyls
        if len(hydroxyl_matches) < (c_count - 1):
            return False, f"Open-chain form lacks sufficient hydroxyl groups, found {len(hydroxyl_matches)} hydroxyls"
    else:
        # Ensure ring carbons have hydroxyl groups
        ring_hydroxyl_count = 0
        ring_atoms = set()
        if furanose_matches:
            ring_atoms.update(furanose_matches[0])
        elif pyranose_matches:
            ring_atoms.update(pyranose_matches[0])
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
                        ring_hydroxyl_count += 1
        if ring_hydroxyl_count < (len(ring_atoms) - 1):
            return False, f"Cyclic form lacks sufficient hydroxyl groups on ring carbons, found {ring_hydroxyl_count} hydroxyls"

    return True, "Molecule is a monosaccharide based on structure"