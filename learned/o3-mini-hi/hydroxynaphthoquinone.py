"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety (a naphthalene-based core carrying two carbonyl groups) is substituted 
            by at least one hydroxy group.
A valid hydroxynaphthoquinone is defined here as a molecule that has:
1. A fused aromatic ring system formed by two six-membered rings sharing exactly two atoms (i.e., a naphthalene core with 10 atoms).
2. Within that fused system, at least two carbonyl (C=O) groups attached directly to the fused ring atoms.
3. And at least one hydroxy substituent (-OH) attached directly.
Note: This approach is heuristic and may not capture every edge‐case.
"""

from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone.
    
    The function first identifies aromatic six-membered rings, then looks for
    pairs that share exactly two atoms (i.e. a naphthalene system, which should
    yield 10 total ring atoms). Finally, for a candidate naphthalene core the code
    checks for at least two carbonyl groups (a double bond from a ring atom to an O)
    and at least one hydroxy group (an -OH group, where the oxygen carries at least one hydrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a hydroxynaphthoquinone, False otherwise.
        str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    # Select aromatic six-membered rings
    rings6 = []
    for ring in ring_info:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            rings6.append(set(ring))
    if not rings6:
        return False, "No aromatic six-membered rings found"
    
    # Look for candidate naphthalene cores (fused system of two 6-membered rings)
    candidate_cores = []
    n_rings = len(rings6)
    for i in range(n_rings - 1):
        for j in range(i + 1, n_rings):
            # They must share exactly 2 atoms to be properly fused as in naphthalene
            inter = rings6[i].intersection(rings6[j])
            if len(inter) == 2:
                # The union should have exactly 10 atoms if it is a naphthalene system
                core = rings6[i].union(rings6[j])
                if len(core) == 10:
                    candidate_cores.append(core)
    
    if not candidate_cores:
        return False, "No fused aromatic ring system resembling a naphthalene (10-atom core) found"
    
    # For each candidate core, count substituents (carbonyl and hydroxy) attached *directly* to the core atoms.
    for core in candidate_cores:
        carbonyl_count = 0
        hydroxy_count = 0
        # For each atom in the core, check its bonds leading out of the core.
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                # Only evaluate bonds to atoms not in the core.
                if neighbor.GetIdx() in core:
                    continue
                # Check for a carbonyl group: a double bond from a ring atom to an oxygen
                if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    carbonyl_count += 1
                # Check for a hydroxy group: a single bond to an oxygen that carries at least one hydrogen.
                if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 8:
                    # Check if the oxygen atom has any bound hydrogens.
                    if neighbor.GetTotalNumHs() > 0:
                        hydroxy_count += 1
        # If this candidate core meets the criteria, we classify the molecule as a hydroxynaphthoquinone.
        if carbonyl_count >= 2 and hydroxy_count >= 1:
            return True, (f"Found a naphthalene core (10 atoms) with {carbonyl_count} carbonyl group(s) "
                          f"and {hydroxy_count} hydroxy substituent(s) attached")
    
    # If no candidate core meets all criteria.
    return False, ("Fused naphthalene core(s) found but none with both at least 2 carbonyl group(s) and at least 1 hydroxy substituent attached")

# Testing examples (optional):
if __name__ == '__main__':
    test_smiles_list = [
        # A few examples provided by the user –
        "O[C@H](C)(C)C1=CC(=O)c2ccccc2C1=O",  # lawsone-like (should be hydroxynaphthoquinone)
        "Oc1cccc2C(=O)C=CC(=O)c12",            # juglone-like (should be hydroxynaphthoquinone)
        "Cc1cc(O)c2C(=O)C=CC(=O)c2c1",          # Ramentaceone-like (should be hydroxynaphthoquinone)
        "COC1=C(C)C(=O)c2c(O)cc(OC\\C=C(/C)CCC=C(C)C)cc2C1=O",  # 7-O-geranyl-2-O,3-dimethylflaviolin (likely not a simple naphthoquinone core)
    ]
    for smi in test_smiles_list:
        result, reason = is_hydroxynaphthoquinone(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")