"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety is substituted 
            by at least one hydroxy group.
A valid hydroxynaphthoquinone is defined here as a molecule that has:
1. A fused aromatic ring system formed by two six-membered rings (a naphthalene-like system)
2. Within that fused system, there are at least two carbonyl (C=O) groups attached directly 
   (as substituents on ring atoms)
3. And at least one hydroxy group (-OH) attached directly.
Note: This approach is heuristic and may not capture every edge‐case.
"""

from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone.

    This function first identifies all aromatic six-membered rings. It then looks
    for pairs of rings that share at least two atoms (indicating a fused system, 
    similar to naphthalene). For each fused ring system, the code examines substituents 
    attached directly to its atoms. A carbonyl group is defined as a double bond to an 
    oxygen atom and a hydroxy substituent is defined as a single bonded oxygen that 
    carries at least one hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a hydroxynaphthoquinone, False otherwise.
        str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    # Filter rings: We require six-membered aromatic rings.
    rings6 = []
    for ring in ring_info:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            rings6.append(set(ring))
    if not rings6:
        return False, "No aromatic six-membered rings found"

    # Look for fused rings (two rings that share at least 2 atoms)
    fused_systems = []
    for i in range(len(rings6) - 1):
        for j in range(i+1, len(rings6)):
            if len(rings6[i].intersection(rings6[j])) >= 2:
                fused = rings6[i].union(rings6[j])
                fused_systems.append(fused)
    if not fused_systems:
        return False, "No fused aromatic ring system (naphthalene-like moiety) found"
    
    # For each fused ring system candidate, count substituents attached
    # outside the fused ring atoms.
    for fused in fused_systems:
        carbonyl_count = 0
        hydroxy_count = 0
        # Check each atom in the fused ring system.
        for idx in fused:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                # Only check substituents not part of the fused ring core.
                if neighbor.GetIdx() in fused:
                    continue

                # Check for carbonyl group: a double bond to oxygen.
                if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    # Additional check: Often in carbonyls the oxygen is terminal.
                    carbonyl_count += 1

                # Check for hydroxy group: single bond to an oxygen with at least one hydrogen.
                if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 8:
                    if neighbor.GetTotalNumHs() > 0:
                        hydroxy_count += 1
        # If criteria are met in this fused system, classify as hydroxynaphthoquinone.
        if carbonyl_count >= 2:
            if hydroxy_count >= 1:
                return True, (f"Found a fused aromatic ring system with {carbonyl_count} carbonyl group(s) "
                              f"and {hydroxy_count} hydroxy substituent(s)")
            else:
                # If one system lacks the OH, we continue checking the other fused systems.
                continue

    # If none of the fused systems meet both criteria.
    return False, ("Fused aromatic ring system(s) found but none with both at least 2 carbonyl "
                   "group(s) and at least 1 hydroxy substituent attached")

# For testing purposes:
if __name__ == '__main__':
    # Example SMILES from the user's list
    test_smiles_list = [
        "O[C@@H](C)(C)C1=CC(=O)c2ccccc2C1=O",  # lawsone-like
        "Oc1cccc2C(=O)C=CC(=O)c12",             # juglone-like
        "Oc1ccc(O)c2C(=O)C=CC(=O)c12",           # naphthazarin-like
        "Cc1cc(O)c2C(=O)C=CC(=O)c2c1",           # Ramentaceone-like
        "COC1=C(C)C(=O)c2c(O)cc(OC\\C=C(/C)CCC=C(C)C)cc2C1=O",  # 7-O-geranyl-2-O,3-dimethylflaviolin-like
    ]
    for smi in test_smiles_list:
        result, reason = is_hydroxynaphthoquinone(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")