"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 4'-hydroxyflavanones
Definition: “Any hydroxyflavanone having a hydroxy substituent located at position 4'.”
Revised approach:
  (1) Identify the flavanone core (chroman-4-one scaffold) via a SMARTS pattern.
  (2) Identify the candidate B ring: a six-membered aromatic ring not fused to the core.
  (3) Verify that the candidate ring is attached by exactly one single bond from an sp³ (non‐aromatic) atom in the core.
  (4) Within the candidate ring, use a graph search to locate the atom at the para position (three bonds away from the attachment)
      and check that it bears an –OH substituent.
Any molecule passing all tests is classified as a 4'-hydroxyflavanone.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if the molecule represented by `smiles` is a 4'-hydroxyflavanone.
    
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
    # This SMARTS is chosen to capture the typical fused benzopyranone system:
    core_smarts = Chem.MolFromSmarts("O=C1CCOc2ccccc21")
    core_matches = mol.GetSubstructMatches(core_smarts)
    if not core_matches:
        return False, ("No appropriate flavanone core was detected – the molecule lacks the fused benzopyranone (chroman-4-one) scaffold typical of flavanones.")
    # Use the first match as representative of the core.
    core_atoms = set(core_matches[0])
    
    # STEP 2: Identify candidate B rings.
    # The candidate B ring is expected to be a six-membered ring that is:
    #    - Fully aromatic
    #    - Not fused with the core (i.e. has no atoms in common with the core)
    #    - Linked to the core by exactly one bond from an sp3 (non‐aromatic) atom.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings detected in the molecule."
    
    candidate_B_rings = []
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        # Exclude rings that share atoms with the core.
        if any(idx in core_atoms for idx in ring):
            continue
        # Enforce full aromaticity.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Now check the connectivity between this ring and the core.
        connecting_bonds = 0
        core_attachment_atom = None   # Atom in the core connected to the B ring.
        b_ring_attachment_idx = None  # Atom in the B ring that connects to core.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    # Ensure that the bond linking the rings is a SINGLE bond.
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue
                    connecting_bonds += 1
                    core_attachment_atom = nbr
                    b_ring_attachment_idx = idx
        # Expect exactly one connecting bond.
        if connecting_bonds != 1:
            continue
        # The attachment from the core should be an sp3 (non‐aromatic) carbon.
        if core_attachment_atom is None or core_attachment_atom.GetHybridization().name != "SP3":
            continue
        candidate_B_rings.append( (tuple(ring), b_ring_attachment_idx) )
    if not candidate_B_rings:
        return False, ("No appropriate six-membered aromatic B ring was detected that is attached by a single bond from an sp³ atom of the flavanone core.")
    
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

    # Helper: Compute shortest path distance within the ring graph between two atoms (BFS).
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
    # In a benzene ring the para position is defined as the atom 3 bonds away from the attachment atom.
    for (ring, attach_idx) in candidate_B_rings:
        g = build_ring_graph(ring)
        found_para_OH = False
        for idx in ring:
            d = ring_distance(g, attach_idx, idx)
            if d != 3:
                continue
            # Check neighbors off the ring atom idx for an -OH group.
            for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # Only consider substituents off the ring.
                # The substituent must be an oxygen with at least one hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    found_para_OH = True
                    break
            if found_para_OH:
                break
        if found_para_OH:
            return True, ("Molecule has a flavanone core (chroman-4-one scaffold) and a non-fused aromatic B ring attached via a single bond from an sp³ atom – the B ring carries an –OH group at the para (4') position relative to its attachment.")
    
    return False, ("The molecule has a flavanone core and a candidate B ring attached appropriately, but no –OH substituent was detected at the para (4') position on the B ring.")

# (Optional) Testing examples.
if __name__ == '__main__':
    # A few test-case examples are provided (names are for clarity).
    test_examples = [
        ("(+)-taxifolin", "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1"),
        ("(2S)-2'-methoxykurarinone", "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1"),
        ("butin", "Oc1ccc2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1"),
        ("naringenin", "Oc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"),
        ("Non-flavanone (acetophenone)", "CC(=O)c1ccccc1")
    ]
    for name, smi in test_examples:
        result, reason = is_4__hydroxyflavanones(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")