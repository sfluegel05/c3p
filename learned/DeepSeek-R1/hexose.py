"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import Mol, MolFromSmiles
from rdkit.Chem.rdchem import AtomValenceException, KekulizeException

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group (aldohexose)
    or a ketone group (ketohexose) in its linear form.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES"

        # Check for exactly six carbons
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count != 6:
            return False, f"Expected 6 carbons, found {c_count}"

        # Check for glycosidic bonds (ethers not in rings)
        # Oxygen connected to two carbons and not in a ring
        ether_pattern = Chem.MolFromSmarts('[OX2;!R]([#6])[#6]')
        if mol.HasSubstructMatch(ether_pattern):
            return False, "Glycosidic bond present"

        # Check for aldehyde group (CH=O)
        aldehyde_pattern = Chem.MolFromSmarts('[CH]=O')
        if mol.HasSubstructMatch(aldehyde_pattern):
            return True, "Aldehyde group detected"

        # Check for ketone group (C=O not in acid/amide)
        ketone_pattern = Chem.MolFromSmarts('[CX3](=O)[#6]')
        ketone_matches = mol.GetSubstructMatches(ketone_pattern)
        # Exclude carboxylic acids and amides
        valid_ketone = False
        for match in ketone_matches:
            atom = mol.GetAtomWithIdx(match[0])
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            # Check if adjacent to O or N (possible acid/amide)
            if 8 not in neighbors and 7 not in neighbors:
                valid_ketone = True
                break
        if valid_ketone:
            return True, "Ketone group detected"

        # Check for cyclic structure with ring oxygen and hydroxyls
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        has_ring_with_oxygen = False
        for ring in rings:
            for atom_idx in ring:
                if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 8:
                    has_ring_with_oxygen = True
                    break
            if has_ring_with_oxygen:
                break
        if has_ring_with_oxygen:
            # Count hydroxyl groups (O with at least one H)
            oh_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0)
            if oh_count >= 4:  # Typical for hexoses
                return True, "Cyclic form with ring oxygen and hydroxyls"
            else:
                return False, f"Cyclic but only {oh_count} hydroxyls"
        else:
            return False, "No aldehyde, ketone, or cyclic structure with ring oxygen"

    except (AtomValenceException, KekulizeException, ValueError):
        return False, "Invalid structure"