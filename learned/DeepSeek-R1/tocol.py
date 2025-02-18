"""
Classifies: CHEBI:39437 tocol
"""
"""
Classifies: CHEBI:46649 tocol
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, GetSSSR

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    Tocols have a chroman-6-ol skeleton with a 13-carbon substituent (saturated or triply unsaturated).
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Ensure ring info is computed
    GetSSSR(mol)
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()

    # Find chroman-6-ol core (fused benzene + dihydropyran with OH)
    chromanol_found = False
    hydroxyl_pos = -1
    for ring in rings:
        if len(ring) == 6:  # Check for 6-membered rings
            oxygen_count = sum(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
            if oxygen_count == 1:  # Dihydropyran ring
                # Check adjacent benzene ring
                for other_ring in rings:
                    if len(other_ring) == 6 and len(set(ring).intersection(other_ring)) >= 2:
                        # Verify hydroxyl group on benzene
                        for atom in other_ring:
                            if (mol.GetAtomWithIdx(atom).GetAtomicNum() == 8 and
                                mol.GetAtomWithIdx(atom).GetTotalNumHs() >= 1):
                                hydroxyl_pos = atom
                                chromanol_found = True
                                break
    if not chromanol_found:
        return False, "No chroman-6-ol core"

    # Find substituent at position 2 (adjacent to oxygen in dihydropyran)
    substituent = None
    for atom in mol.GetAtomWithIdx(hydroxyl_pos).GetNeighbors():
        if atom.GetAtomicNum() == 6:  # Carbon adjacent to hydroxyl
            substituent = atom
            break
    if not substituent:
        return False, "No substituent at position 2"

    # Traverse the substituent chain
    chain = []
    visited = set()
    stack = [(substituent, 0)]  # (atom, bond count)
    while stack:
        atom, bonds = stack.pop()
        if atom.GetIdx() in visited:
            continue
        visited.add(atom.GetIdx())
        chain.append(atom)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    stack.append((neighbor, bonds + 1))

    # Check chain properties: length and unsaturation
    carbon_count = len(chain)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE and
                       bond.GetBeginAtom().GetIdx() in [a.GetIdx() for a in chain] and
                       bond.GetEndAtom().GetIdx() in [a.GetIdx() for a in chain])
    
    if carbon_count >= 13 and (double_bonds == 0 or double_bonds == 3):
        return True, "Chroman-6-ol core with 13+ carbon chain (0 or 3 double bonds)"
    return False, f"Chain has {carbon_count} carbons and {double_bonds} double bonds (needs 13+ and 0/3)"