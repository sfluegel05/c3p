"""
Classifies: CHEBI:46640 diketone
"""
"""
Classifies: A diketone – a compound that contains exactly two independent ketone functionalities.
A “ketone functionality” (as used here) is defined as a carbonyl group R–CO–R,
where the carbonyl carbon:
  • Is sp2‐hybridized,
  • Has exactly three neighbors (one double-bonded oxygen and two single-bonded carbons),
  • Has no attached hydrogens,
and if the two candidate ketone carbons lie in a common ring that is mostly aromatic
(quinone‐like), then they are not counted.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_diketone(smiles: str):
    """
    Determines if a molecule is a diketone (contains exactly two independent ketone functionalities)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if exactly two independent ketone functionalities are found (and they are not both in a quinone-like aromatic ring),
              False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule (this computes ring info, hybridizations, etc.)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Could not sanitize molecule: " + str(e)
    
    rings = mol.GetRingInfo().AtomRings()  # each ring is a tuple of atom indices

    # Helper function to check if an atom lies in any ring that is fully aromatic.
    # (Used later when deciding if two candidate ketone carbons share a quinone-like system.)
    def in_fully_aromatic_ring(atom_idx):
        for ring in rings:
            if atom_idx in ring:
                # Determine aromatic fraction of the ring:
                aromatic_atoms = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
                if len(ring) >= 5 and (aromatic_atoms / len(ring)) > 0.5:
                    return True
        return False

    # Identify candidate ketone carbon atoms using our revised criteria.
    candidate_atoms = []
    for atom in mol.GetAtoms():
        # We are interested only in carbon atoms.
        if atom.GetSymbol() != "C":
            continue
        # For ketone carbon, we expect exactly 3 neighbors.
        if atom.GetDegree() != 3:
            continue
        # Must have no attached hydrogens.
        if atom.GetTotalNumHs() != 0:
            continue

        # We count bonds: exactly one double bond (to an oxygen) and two single bonds (to carbons)
        doubleO_count = 0
        singleC_count = 0
        valid = True  # flag to track if any bond fails our criteria
        for bond in atom.GetBonds():
            btype = bond.GetBondType()
            nbr = bond.GetOtherAtom(atom)
            if btype == rdchem.BondType.DOUBLE:
                if nbr.GetAtomicNum() == 8:
                    doubleO_count += 1
                else:
                    valid = False
                    break
            elif btype == rdchem.BondType.SINGLE:
                # Only count if neighbor is carbon.
                if nbr.GetAtomicNum() == 6:
                    singleC_count += 1
                else:
                    valid = False
                    break
            else:
                valid = False
                break
        if not valid:
            continue

        if doubleO_count == 1 and singleC_count == 2:
            candidate_atoms.append(atom)
    
    ketone_count = len(candidate_atoms)
    if ketone_count != 2:
        return False, f"Found {ketone_count} ketone functionalities; exactly 2 are required for a diketone"
    
    # Apply the quinone-like filter.
    # If both candidate ketone carbons lie in a common ring that is mostly aromatic, we assume a quinone-like system.
    candidate_idxs = {atom.GetIdx() for atom in candidate_atoms}
    for ring in rings:
        if candidate_idxs.issubset(ring):
            aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
            if len(ring) >= 5 and (aromatic_count / len(ring)) > 0.5:
                return False, "Ketone functionalities appear as part of a quinone-like aromatic system"
    
    return True, "Contains exactly two ketone functionalities"

# Example usage (optional test cases):
if __name__ == "__main__":
    test_cases = [
        ("O=C(CCCCCCCC)C(=O)C", "2,3-Undecanedione"),
        ("O=C(C(CCCCCCCC)CC)CC(=O)CCC", "7-Ethylpentadecane-4,6-dione"),
        ("O=C(CCCCCCC)C(=O)C", "2,3-Decanedione"),
        ("CCCC(C(C)=O)=O", "2,3-Heptanedione"),
        ("[H][C@@]1(CC[C@@]2(C)C3=C(CC[C@]12C)[C@@]1(C)CCC(=O)C(C)(C)[C@]1([H])CC3=O)[C@H](C)CC/C=C(/C)CO", "ganoderone A"),
        ("OC1ccc(cc1)C(=O)CC(=O)c1ccc(O)cc1O", "licodione"),
        # Additional examples can be added here.
    ]
    
    for smi, name in test_cases:
        res, reason = is_diketone(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {res}\nReason: {reason}\n{'-'*60}")