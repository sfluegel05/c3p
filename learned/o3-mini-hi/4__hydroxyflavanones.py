"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 4'-hydroxyflavanones
Definition: “Any hydroxyflavanone having a hydroxy substituent located at position 4'.”
Revised approach:
  (1) Identify the flavanone core as a chroman-4-one. The SMARTS pattern captures the fused ring system 
      (the benzopyranone core) of a typical flavanone.
  (2) Identify a candidate B ring: a six-membered fully aromatic ring that is not fused to the core but attached 
      via a single bond from an sp³ (non‐aromatic) core atom.
  (3) In the candidate B ring, verify that the atom directly attached to the core has an –OH substituent at the para position.
Any molecule passing all tests is classified as a 4'-hydroxyflavanone.
"""
from rdkit import Chem
from collections import deque

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if the molecule represented by `smiles` is a 4'-hydroxyflavanone.
    The approach is to detect a flavanone (chroman-4-one) core and then inspect the appended B ring.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a 4'-hydroxyflavanone, False otherwise.
       str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Identify the flavanone (chroman-4-one) core.
    # The following SMARTS is intended to capture the classical 2-phenylchroman-4-one scaffold.
    # It looks for a 6-membered ring incorporating a ketone (C=O) and an ether O.
    core_smarts = Chem.MolFromSmarts("O=C1CCOc2ccccc21")
    core_matches = mol.GetSubstructMatches(core_smarts)
    if not core_matches:
        return False, ("No appropriate flavanone core was detected – the molecule lacks the fused benzopyranone (chroman-4-one) scaffold typical of flavanones.")
    
    # For further analysis, use the first matched core.
    core_atoms = set(core_matches[0])
    
    # STEP 2: Identify candidate B rings.
    # In flavanones, the B ring is typically the phenyl substituent attached at C2 of the chromanone.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings detected in the molecule."
    
    candidate_B_rings = []
    for ring in atom_rings:
        # We are interested in six-membered rings.
        if len(ring) != 6:
            continue
        # Exclude rings that are fused with (i.e. share atoms with) the core.
        if any(idx in core_atoms for idx in ring):
            continue
        # The ring must be fully aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Count bonds connecting this ring to the core.
        connecting_bonds = 0
        core_attachment_atom = None   # atom in core that connects to B ring
        b_ring_attachment_idx = None  # atom in B ring that connects to core
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    connecting_bonds += 1
                    core_attachment_atom = nbr
                    b_ring_attachment_idx = idx
        # Expect exactly one connection.
        if connecting_bonds != 1:
            continue
        # In a typical flavanone the core atom contacting the B ring is sp3 (non‐aromatic).
        if core_attachment_atom is None or core_attachment_atom.GetHybridization().name != "SP3":
            continue
        candidate_B_rings.append( (tuple(ring), b_ring_attachment_idx) )
    if not candidate_B_rings:
        return False, ("No appropriate six-membered aromatic B ring was detected that is attached via a single bond "
                       "from an sp³ atom in the flavanone core.")
    
    # Helper: Build an undirected graph (adjacency list) for atoms in a given ring.
    def build_ring_graph(ring):
        g = {i: set() for i in ring}
        ring_set = set(ring)
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                g[a1].add(a2)
                g[a2].add(a1)
        return g

    # Helper: Compute shortest path distance within the ring graph between two atoms (using BFS).
    def ring_distance(g, start, end):
        visited = {start}
        queue = deque([(start, 0)])
        while queue:
            current, d = queue.popleft()
            if current == end:
                return d
            for nbr in g[current]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, d+1))
        return None

    # STEP 3: For each candidate B ring, check for an –OH group at the para position relative to the attachment.
    # In a benzene ring, the para position is 3 bonds away.
    for (ring, attach_idx) in candidate_B_rings:
        g = build_ring_graph(ring)
        # Within the B ring, the attachment atom is already defined (attach_idx).
        # Search for an atom in the ring with an external hydroxyl (-OH) that sits exactly 3 bonds away (para).
        found_para_OH = False
        for idx in ring:
            # Check ring distance from the attachment atom.
            d = ring_distance(g, attach_idx, idx)
            if d != 3:
                continue
            # Look at substituents off the ring atom.
            for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if the neighbor is oxygen and if it carries at least one hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    found_para_OH = True
                    break
            if found_para_OH:
                break
        if found_para_OH:
            return True, ("Molecule has a flavanone core (chroman-4-one scaffold) and a non-fused aromatic B ring, "
                          "where the B ring carries an –OH at the para (4') position relative to its attachment.")
    
    return False, ("The molecule contains a flavanone core and a candidate B ring, "
                   "but no –OH substituent was found at the para (4') position on the B ring.")

# (Optional) Testing examples.
if __name__ == '__main__':
    test_examples = [
        ("(+)-taxifolin", "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1"),
        ("(2S)-2'-methoxykurarinone", "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1"),
        ("butin", "Oc1ccc2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1"),
        ("naringenin", "Oc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"),
        ("non-flavanone (acetophenone)", "CC(=O)c1ccccc1")
    ]
    for name, smi in test_examples:
        result, reason = is_4__hydroxyflavanones(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")