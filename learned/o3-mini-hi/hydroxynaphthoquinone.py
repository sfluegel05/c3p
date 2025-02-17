"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthoquinone moiety 
            (a fused bicyclic ring system that is "naphthalene‐like" comprising 10 carbon atoms,
             two of which bear a double‐bonded oxygen) is substituted by at least one hydroxy group.
            
A valid hydroxynaphthoquinone is defined here as having:
1. A fused bicyclic ring system formed by two 6–membered rings sharing exactly 2 atoms 
   (i.e. a candidate naphthalene core consisting of 10 atoms).
2. Every atom in the candidate core must be carbon (atomic number 6).
3. At least two carbonyl groups (a double bond from a core carbon to an oxygen outside the core).
4. And at least one hydroxy substituent (-OH group) attached directly to a core carbon.
Note: This method is heuristic and may miss some edge–cases.
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone.
    
    The function:
    1. Parses the SMILES string.
    2. Uses ring information to find candidate fused bicyclic (naphthalene-like) systems.
    3. For each candidate of exactly 10 atoms, checks that all atoms are carbon.
    4. For each candidate, counts carbonyl groups (direct double bonds to an oxygen, not in the core)
       and hydroxy substituents (single bonds to oxygen that carries at least one attached hydrogen).
    5. If a candidate has at least two such carbonyl groups and at least one hydroxy, it is classified as a hydroxynaphthoquinone.
    
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
    rings6 = []
    # Collect all rings having 6 atoms (do not force aromatic flag, because carbonyl groups sometimes unset aromaticity)
    for ring in ring_info:
        if len(ring) == 6:
            rings6.append(set(ring))
    if not rings6:
        return False, "No six‐membered rings found in molecule"
    
    candidate_cores = []
    n_rings = len(rings6)
    # Look for pairs of rings that share exactly 2 atoms and whose union forms a 10–atom fused core.
    for i in range(n_rings - 1):
        for j in range(i + 1, n_rings):
            shared = rings6[i].intersection(rings6[j])
            if len(shared) == 2:
                core = rings6[i].union(rings6[j])
                if len(core) == 10:
                    # Check if every atom in the core is carbon (atomic number 6)
                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in core):
                        candidate_cores.append(core)
                        
    if not candidate_cores:
        return False, "No fused 10–atom (naphthalene-like) carbon core found"

    # Now evaluate each candidate core
    for core in candidate_cores:
        carbonyl_count = 0
        hydroxy_count = 0
        # For each atom in the core:
        for idx in core:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                # Only consider bonds going out of the core
                if neighbor.GetIdx() in core:
                    continue
                # Check for a carbonyl group:
                # It should be a double bond from a core carbon to an oxygen (outside the core).
                if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
                    carbonyl_count += 1
                # Check for a hydroxy group:
                # A single bond from a core carbon to oxygen that carries at least one hydrogen.
                if bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetAtomicNum() == 8:
                    # Sometimes the number of implicit/explicit hydrogens is not updated.
                    # We ask for total hydrogens on the oxygen.
                    if neighbor.GetTotalNumHs() > 0:
                        hydroxy_count += 1
        
        # Now require at least two carbonyl groups and at least one hydroxy group
        if carbonyl_count >= 2 and hydroxy_count >= 1:
            return True, (f"Found a naphthalene-like 10–atom carbon core with {carbonyl_count} carbonyl group(s) "
                          f"and {hydroxy_count} hydroxy substituent(s) attached")
    
    return False, ("Fused naphthalene-like core(s) were found, but none with both at least 2 carbonyl group(s) "
                   "and at least 1 hydroxy substituent attached")

# (Optional) Testing examples – uncomment the following lines to run some tests:
# if __name__ == '__main__':
#     test_smiles_list = [
#         "O[C@H](C)(C)C1=CC(=O)c2ccccc2C1=O",  # lawsone-like: expected True
#         "Oc1cccc2C(=O)C=CC(=O)c12",            # juglone-like: expected True
#         "Cc1cc(O)c2C(=O)C=CC(=O)c2c1",          # Ramentaceone-like: expected True
#         "COC1=C(C)C(=O)c2c(O)cc(OC\\C=C(/C)CCC=C(C)C)cc2C1=O",  # 7-O-geranyl-2-O,3-dimethylflaviolin: may be rejected by our filter
#     ]
#     for smi in test_smiles_list:
#         result, reason = is_hydroxynaphthoquinone(smi)
#         print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")