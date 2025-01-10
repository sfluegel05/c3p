"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as any member of class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage to a carbohydrate part (usually a mono-, di- or tri-saccharide).
    Some bacterial glycolipids have the sugar part acylated by one or more fatty acids and the glycerol part may be absent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to check if molecule contains a sugar moiety
    found_sugar = False

    # Detect sugar rings (5 or 6-membered rings with oxygen, carbons having hydroxyl groups)
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        ring_size = len(ring)
        if ring_size == 5 or ring_size == 6:
            oxygens_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            # Consider rings with 1 oxygen in ring
            if oxygens_in_ring == 1:
                carbons_with_oh = 0
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 6:  # Carbon atom
                        # Check if carbon is sp3 hybridized
                        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                            break
                        # Check for attached hydroxyl group
                        hydroxyl = False
                        for nbr in atom.GetNeighbors():
                            if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                                hydroxyl = True
                        if hydroxyl:
                            carbons_with_oh +=1
                        else:
                            break
                else:
                    # All carbons in ring have hydroxyl groups
                    found_sugar = True
                    break  # Found a sugar ring

    if not found_sugar:
        return False, "No carbohydrate moiety detected"

    # Check for long aliphatic chains (lipid moiety)
    lipid_found = False
    num_carbons_in_chains = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1:
            chain_length = 1
            visited = set()
            stack = [(atom, 0)]
            while stack:
                current_atom, depth = stack.pop()
                if current_atom.GetIdx() in visited:
                    continue
                visited.add(current_atom.GetIdx())
                if current_atom.GetAtomicNum() == 6:
                    chain_length += 1
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() in [6, 1]:  # Carbon or hydrogen
                        stack.append((neighbor, depth + 1))
            if chain_length >= 8:
                lipid_found = True
                break

    if not lipid_found:
        return False, "No lipid moiety detected"

    # Check for glycosidic linkage between sugar and lipid
    # Look for C-O-C linkage between sugar ring and lipid chain
    glycosidic_linkage_found = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6 or atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8:
                # Check if one atom is in sugar ring and the other is outside
                idx1 = atom1.GetIdx()
                idx2 = atom2.GetIdx()
                in_sugar_ring1 = any(idx1 in ring for ring in rings)
                in_sugar_ring2 = any(idx2 in ring for ring in rings)
                if in_sugar_ring1 != in_sugar_ring2:
                    glycosidic_linkage_found = True
                    break

    if not glycosidic_linkage_found:
        return False, "No glycosidic linkage between sugar and lipid moiety found"

    return True, "Contains carbohydrate moiety linked to lipid via glycosidic bond"