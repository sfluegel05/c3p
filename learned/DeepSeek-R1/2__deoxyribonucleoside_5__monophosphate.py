"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate.
    Criteria: Deoxyribose (no 2'-OH), phosphate at 5' position, nucleobase at 1' position.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # General SMARTS for deoxyribose core with phosphate at C5' and base at C1'
    # More flexible stereochemistry
    pattern = Chem.MolFromSmarts(
        "[O;R1]"
        "C(-[C@H](O)C(-COP(=O)([O-])[O-])[C@H]1O)"  # C1' connected to ring O, C2', C3', C4', C5'
        "1"  # Ring closure
        "[N,n,o,O]"  # Base attached to C1'
    )
    if mol.HasSubstructMatch(pattern):
        # Verify C2' has no -OH
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            c2_prime = None
            # Get C1' atom index from match (index 1 in SMARTS)
            c1_prime_idx = match[1]
            c1_prime = mol.GetAtomWithIdx(c1_prime_idx)
            # Find C2' (connected to C1' in the ring)
            for neighbor in c1_prime.GetNeighbors():
                if neighbor.GetIdx() in match:  # Check if in the ring part
                    bond = mol.GetBondBetweenAtoms(c1_prime_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.BondType.SINGLE:
                        c2_prime = neighbor
                        break
            if c2_prime:
                # Check for absence of -OH on C2'
                has_oh = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0 
                            for atom in c2_prime.GetNeighbors())
                if not has_oh:
                    return True, "Core structure match with 5'-phosphate and no 2'-OH"
        
    # Manual verification for edge cases
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        # Check for furanose oxygen
        o_in_ring = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
        if not o_in_ring:
            continue

        # Find C1' (connected to base and ring O)
        c1_prime = None
        for a in ring:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetSymbol() != 'C':
                continue
            # Connected to ring O and has a heteroatom (base) outside the ring
            neighbors = atom.GetNeighbors()
            o_neighbor = any(n.GetAtomicNum() == 8 and n.GetIdx() in ring for n in neighbors)
            base_neighbor = any(n.GetAtomicNum() in [7,8] and n.GetIdx() not in ring for n in neighbors)
            if o_neighbor and base_neighbor:
                c1_prime = atom
                break
        if not c1_prime:
            continue

        # Check C2' (next to C1' in ring) has no OH
        c2_prime = None
        for a in ring:
            if a == c1_prime.GetIdx():
                continue
            if mol.GetBondBetweenAtoms(c1_prime.GetIdx(), a):
                c2_prime = mol.GetAtomWithIdx(a)
                break
        if not c2_prime or any(n.GetAtomicNum() == 8 and n.GetTotalNumHs() > 0 for n in c2_prime.GetNeighbors()):
            continue

        # Check C5' (connected to C4') has phosphate
        c4_prime = None
        for a in ring:
            if a == c1_prime.GetIdx() or a == c2_prime.GetIdx():
                continue
            if mol.GetBondBetweenAtoms(c1_prime.GetIdx(), a):
                # Assuming ring order, find C4'
                c4_prime = mol.GetAtomWithIdx(a)
                break
        if not c4_prime:
            continue

        c5_prime = None
        for neighbor in c4_prime.GetNeighbors():
            if neighbor.GetIdx() not in ring and neighbor.GetSymbol() == 'C':
                c5_prime = neighbor
                break
        if not c5_prime:
            continue

        # Check phosphate attached to C5'
        phosphate = False
        for bond in c5_prime.GetBonds():
            other = bond.GetOtherAtom(c5_prime)
            if other.GetSymbol() == 'O':
                for o_bond in other.GetBonds():
                    if o_bond.GetOtherAtom(other).GetSymbol() == 'P':
                        phosphate = True
                        break
                if phosphate:
                    break
        if phosphate:
            return True, "Manual check: deoxyribose, 5'-phosphate, base"

    return False, "Does not meet criteria"