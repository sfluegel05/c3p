"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import Mol, MolFromSmiles, rdMolDescriptors
from rdkit.Chem.rdchem import AtomValenceException

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

        # Basic checks for monosaccharide structure
        # Total carbons must be exactly 6
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        if c_count != 6:
            return False, f"Has {c_count} carbons (needs 6)"

        # Exclude derivatives with non-sugar components
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() not in {1, 6, 8}:  # Only allow H, C, O
                return False, "Contains non-C/H/O atoms"
            if atom.GetFormalCharge() != 0:
                return False, "Charged groups present"

        # Check for forbidden functional groups
        forbidden_patterns = [
            Chem.MolFromSmarts(patt) for patt in [
                '[#15]', '[#16]',          # Phosphorus/Sulfur
                '[NX3][CX3](=O)',          # Amides
                '[OX2][CX3]=O',            # Esters
                '[NX3+!H0]',               # Amines (charged)
                '[CX3](=O)[OX2H1]'         # Carboxylic acids
            ]
        ]
        for patt in forbidden_patterns:
            if mol.HasSubstructMatch(patt):
                return False, "Contains forbidden group"

        # Check for glycosidic bonds (ethers not in rings)
        ether_pattern = Chem.MolFromSmarts('[OX2;!R]([#6])[#6]')
        if mol.HasSubstructMatch(ether_pattern):
            return False, "Glycosidic bond present"

        # Check linear forms first
        # Aldehyde check (terminal CH=O, allowing explicit H)
        aldehyde_pattern = Chem.MolFromSmarts('[CH0X3]=O')  # Matches [H]C=O
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
        if aldehyde_matches:
            for match in aldehyde_matches:
                ald_atom = mol.GetAtomWithIdx(match[0])
                # Must be terminal (only one neighbor)
                if ald_atom.GetDegree() == 1:
                    return True, "Aldehyde group detected (aldohexose)"
            return False, "Aldehyde not terminal"

        # Ketone check (C=O at position 2 in 6-carbon chain)
        ketone_pattern = Chem.MolFromSmarts('[CX3]=O')
        ketone_matches = mol.GetSubstructMatches(ketone_pattern)
        if ketone_matches:
            # Get longest chain (assumed to be the sugar backbone)
            main_chain = Chem.GetLongestChain(mol)
            if len(main_chain) == 6:  # Must be 6 carbons in chain
                ketone_positions = [mol.GetAtomWithIdx(m[0]).GetIdx() for m in ketone_matches]
                # Position 2 in chain is index 1 (0-based)
                if main_chain[1] in ketone_positions:
                    return True, "Ketone at position 2 (ketohexose)"
            return False, "Ketone not at position 2"

        # Check cyclic forms (pyranose/furanose)
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        valid_cyclic = False
        for ring in rings:
            if len(ring) not in {5,6}:
                continue  # Only 5/6 membered rings
            
            # Check for exactly one oxygen in the ring
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count != 1:
                continue
            
            # Count hydroxyl groups attached to ring atoms
            oh_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    continue  # Skip the ring oxygen
                # Check for -OH groups
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() > 0:
                        oh_count += 1
                        break  # Count once per ring atom
            
            # Pyranose (6-membered) needs >=4 OHs, furanose (5) >=3
            if (len(ring) == 6 and oh_count >=4) or (len(ring) ==5 and oh_count >=3):
                valid_cyclic = True
                break
        
        if valid_cyclic:
            return True, "Valid cyclic hexose form"
        
        return False, "No aldehyde/ketone or valid cyclic structure"

    except AtomValenceException:
        return False, "Invalid valence"
    except Exception as e:
        return False, f"Processing error: {str(e)}"