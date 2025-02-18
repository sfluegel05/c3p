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
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    Criteria: Deoxyribose (no 2'-OH), phosphate at 5' position, nucleobase at 1' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # General SMARTS pattern for deoxyribose with 5'-phosphate and any base
    # Allows variable stereochemistry and base attachment
    pattern = Chem.MolFromSmarts(
        "[O;R1]"
        "[C@H]1"  # C1'
        "[C@@H](O)"  # C3' (with hydroxyl)
        "[C@@H](COP(=O)([OH,O-])[OH,O-])"  # C4' connected to C5' with phosphate
        "O"  # O connecting to C1' (part of ring)
        "[C@H]1"  # Completes the ring
        "[N,n,O,o]"  # Base attachment (any heteroatom)
    )
    if mol.HasSubstructMatch(pattern):
        return True, "Core structure with deoxyribose, 5'-phosphate, and base"

    # Alternative pattern with different stereochemistry
    alt_pattern = Chem.MolFromSmarts(
        "[O;R1]"
        "[C@@H]1"  # C1'
        "[C@H](O)"  # C3'
        "[C@H](COP(=O)([OH,O-])[OH,O-])"  # C4'
        "O"  # Ring O
        "[C@@H]1"  # Completes ring
        "[N,n,O,o]"  # Base
    )
    if mol.HasSubstructMatch(alt_pattern):
        return True, "Alternative stereochemistry match"

    # Check for deoxyribose characteristics manually
    # Find all 5-membered rings with oxygen
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 5:
            continue
        o_present = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
        if not o_present:
            continue

        # Check for 2'-deoxy (no OH on C2')
        # Find the ring oxygen and adjacent carbons (C1' and C4')
        for o_idx in [a for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 8]:
            o_atom = mol.GetAtomWithIdx(o_idx)
            neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() in ring]
            if len(neighbors) != 2:
                continue
            c1_prime, c4_prime = neighbors

            # Find C2' as neighbor of C1' in the ring
            c2_prime = None
            for a in ring:
                if a == o_idx or a == c1_prime.GetIdx() or a == c4_prime.GetIdx():
                    continue
                if mol.GetBondBetweenAtoms(c1_prime.GetIdx(), a):
                    c2_prime = mol.GetAtomWithIdx(a)
                    break
            if not c2_prime:
                continue

            # Check C2' has no hydroxyl
            has_oh = any(bond.GetOtherAtom(c2_prime).GetAtomicNum() == 8 and 
                         bond.GetOtherAtom(c2_prime).GetTotalNumHs() > 0 
                         for bond in c2_prime.GetBonds())
            if has_oh:
                continue

            # Check C5' (connected to C4') has phosphate
            c5_prime = None
            for neighbor in c4_prime.GetNeighbors():
                if neighbor.GetIdx() not in ring:
                    c5_prime = neighbor
                    break
            if not c5_prime or c5_prime.GetSymbol() != 'C':
                continue

            # Check C5' is connected to phosphate via O
            phosphate_found = False
            for bond in c5_prime.GetBonds():
                other_atom = bond.GetOtherAtom(c5_prime)
                if other_atom.GetSymbol() == 'O':
                    for o_bond in other_atom.GetBonds():
                        p_atom = o_bond.GetOtherAtom(other_atom)
                        if p_atom.GetSymbol() == 'P':
                            phosphate_found = True
                            break
                if phosphate_found:
                    break
            if not phosphate_found:
                continue

            # Check C1' is connected to a nucleobase (heterocycle)
            base_attached = None
            for neighbor in c1_prime.GetNeighbors():
                if neighbor.GetIdx() not in ring:
                    base_attached = neighbor
                    break
            if not base_attached:
                continue

            # Verify base is a heterocycle (ring with N/O)
            base_rings = ring_info.AtomRings(base_attached.GetIdx())
            if not base_rings:
                continue
            hetero_in_ring = any(
                mol.GetAtomWithIdx(a).GetAtomicNum() in [7, 8] 
                for a in base_rings[0]
            )
            if hetero_in_ring:
                return True, "Manual verification: deoxyribose with 5'-phosphate and nucleobase"

    return False, "Does not meet all criteria"