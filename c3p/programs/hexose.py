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

        # Exclude derivatives with non-sugar components
        # Check for phosphorus (phosphates)
        if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
            return False, "Contains phosphorus"
        # Check for sulfur (sulfates)
        if any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms()):
            return False, "Contains sulfur"
        # Check for amides
        amide_pattern = Chem.MolFromSmarts('[NX3][CX3](=O)')
        if mol.HasSubstructMatch(amide_pattern):
            return False, "Contains amide"
        # Check for esters (O connected to carbonyl)
        ester_pattern = Chem.MolFromSmarts('[OX2][CX3]=O')
        if mol.HasSubstructMatch(ester_pattern):
            return False, "Contains ester"

        # Check for glycosidic bonds (ethers not in rings)
        ether_pattern = Chem.MolFromSmarts('[OX2;!R]([#6])[#6]')
        if mol.HasSubstructMatch(ether_pattern):
            return False, "Glycosidic bond present"

        # Check for aldehyde group (CH=O)
        aldehyde_pattern = Chem.MolFromSmarts('[CH]=O')
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
        if aldehyde_matches:
            # Verify it's terminal in the main chain
            for match in aldehyde_matches:
                ald_atom = mol.GetAtomWithIdx(match[0])
                # Check if aldehyde is at end of chain (only one neighbor)
                if ald_atom.GetDegree() == 1:
                    return True, "Aldehyde group detected"
            return False, "Aldehyde not terminal"

        # Check for ketone group (C=O not in acid/amide)
        ketone_pattern = Chem.MolFromSmarts('[CX3](=O)[#6]')
        ketone_matches = mol.GetSubstructMatches(ketone_pattern)
        valid_ketone = False
        for match in ketone_matches:
            atom = mol.GetAtomWithIdx(match[0])
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if 8 not in neighbors and 7 not in neighbors:  # Not acid/amide
                # Check if ketone is at position 2 in a 6-carbon chain
                # This is simplified check - assumes linear form
                chain = Chem.GetLongestChain(mol)
                if len(chain) >= 6:
                    # Position 2 in chain (0-based index 1)
                    if atom.GetIdx() == chain[1]:
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
            if oh_count >= 3:
                return True, "Cyclic form with ring oxygen and hydroxyls"
            else:
                return False, f"Cyclic but only {oh_count} hydroxyls"
        else:
            return False, "No aldehyde, ketone, or cyclic structure with ring oxygen"

    except (AtomValenceException, KekulizeException, ValueError):
        return False, "Invalid structure"