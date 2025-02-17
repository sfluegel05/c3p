"""
Classifies: CHEBI:33567 catecholamine
"""
#!/usr/bin/env python
"""
Classifies: catecholamine 
Definition: catecholamines include 4-(2-Aminoethyl)pyrocatechol and derivatives.
A catecholamine must contain:
  (i) a catechol moiety—a benzene ring (six-membered aromatic ring) with two adjacent oxygen substituents,
  (ii) an aminoethyl side chain attached directly to the catechol ring (i.e. a two carbon aliphatic chain ending in an amine).
This program uses RDKit for molecule parsing and heuristic substructure detection.
"""

from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    The molecule must have:
      - An aromatic benzene ring (6-membered) that bears at least 2 adjacent oxygen substituents 
        (these are not part of the ring, e.g. –OH, -OCH3, etc.),
      - An aminoethyl side chain (i.e. a 2-carbon chain ending with an amine) attached directly to one of the ring atoms.
    
    Args:
      smiles (str): SMILES string representing the molecule.
    
    Returns:
      bool: True if the molecule is classified as a catecholamine, False otherwise.
      str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    catechol_found = False    # At least one catechol ring found
    chain_attached = False    # At least one catechol ring has an aminoethyl side chain
    msg_no_ring = "No catechol ring (aromatic 6-membered ring with adjacent oxygen substituents) found"
    msg_no_chain = "Catechol ring detected but no attached aminoethyl side chain (2-carbon chain ending in N) found"
    
    # Helper: Given a ring atom index, check if it has an aminoethyl chain attached.
    def has_aminoethyl_chain(ring_atom_idx):
        # For every neighbor outside the ring, try to detect a chain: ring_atom --> carbon (first) --> carbon (second) --> nitrogen
        ring_atom = mol.GetAtomWithIdx(ring_atom_idx)
        for nbr1 in ring_atom.GetNeighbors():
            # Ensure neighbor is not part of the ring.
            if nbr1.GetIdx() == ring_atom_idx:
                continue
            if nbr1.GetIdx() in current_ring: 
                continue
            # First atom should be an aliphatic carbon.
            if nbr1.GetAtomicNum() != 6 or nbr1.GetIsAromatic():
                continue
            # Explore second bond from nbr1
            for nbr2 in nbr1.GetNeighbors():
                if nbr2.GetIdx() == ring_atom_idx:
                    continue
                # We allow the first carbon to be substituted (e.g. -OH attached)
                # Look for a second carbon in the chain (should be aliphatic carbon)
                if nbr2.GetAtomicNum() != 6 or nbr2.GetIsAromatic():
                    continue
                # Now, look for a nitrogen attached to nbr2.
                for nbr3 in nbr2.GetNeighbors():
                    if nbr3.GetIdx() in (nbr1.GetIdx(), ring_atom_idx):
                        continue
                    if nbr3.GetAtomicNum() == 7:
                        return True
        return False

    # Look for catechol ring: an aromatic 6-membered ring with at least two adjacent oxygen substituents.
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Check that all atoms in the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # For each ring atom, check if an oxygen substituent (not in the ring) is bonded via a single bond.
        oxy_substituted = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxy_substituted.add(idx)
                    break
        # To be catechol, we need at least two adjacent ring atoms with oxygen substituents.
        if len(oxy_substituted) < 2:
            continue
        # Check for adjacency in the ring (cyclic order).
        ring_is_catechol = False
        ring_len = len(ring)
        for i in range(ring_len):
            current_atom_idx = ring[i]
            next_atom_idx = ring[(i+1) % ring_len]
            if current_atom_idx in oxy_substituted and next_atom_idx in oxy_substituted:
                ring_is_catechol = True
                break
        if not ring_is_catechol:
            continue
        
        # If this ring qualifies as a catechol moiety, mark that one was found.
        catechol_found = True
        # Save the current ring indices for use in chain search.
        current_ring = set(ring)
        # Now, for each atom in the ring, try to detect the aminoethyl chain attached directly.
        for idx in ring:
            if has_aminoethyl_chain(idx):
                chain_attached = True
                # We break out if we find at least one catechol ring with an attached aminoethyl chain.
                break
        if chain_attached:
            break  # no need to search further

    if not catechol_found:
        return False, msg_no_ring
    if not chain_attached:
        return False, msg_no_chain

    return True, "Contains a catechol moiety with an aminoethyl side chain attached directly to the ring"

# For testing purposes
if __name__ == "__main__":
    # (The examples below are provided for testing and demonstrate that a molecule must have 
    # both the catechol ring and an aminoethyl side chain directly connected)
    test_molecules = [
        # True positives (catecholamine examples)
        ("C(CNCCCCCCNCCC1=CC=CC=C1)C2=CC(O)=C(C=C2)O", "dopexamine"),
        ("OC(=O)C1CC(=C\\C=N/CCc2ccc(O)c(O)c2)/C=C(N1)C(O)=O", "Miraxanthin-V"),
        ("C[C@H](N)[C@H](O)c1ccc(O)c(O)c1", "(-)-alpha-Methylnoradrenaline"),
        ("C=1(C=C(C(O)=CC1)O)CCN.Cl", "Dopamine hydrochloride"),
        ("C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O", "(R)-dobutamine"),
        ("C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1", "(S)-dobutamine"),
        ("OC(=O)[C@H](Cc1ccc(O)c(O)c1)\\N=C/C=C1C[C@H](NC(=C\\1)C(O)=O)C(O)=O", "Dopaxanthin"),
        ("CNC[C@H](O)c1ccc(O)c(O)c1", "(R)-adrenaline"),
        ("CC(N)C(O)c1ccc(O)c(O)c1", "Nordephrine"),
        ("NC[C@@H](O)c1ccc(O)c(O)c1", "(S)-noradrenaline"),
        ("[C@@H]([C@@H](N)C)(O)C1=CC(O)=C(C=C1)O", "alpha-methylnoradrenaline"),
        ("C1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN", "4-(2-aminoethyl)-5-nitrobenzene-1,2-diol"),
        ("C1=CC(=C(C(=C1CCN)O)O)O", "4-(2-aminoethyl)benzene-1,2,3-triol"),
        ("O[C@@H](CNCCCCc1ccc(O)cc1)c1ccc(O)c(O)c1", "arbutamine"),
        ("CNC[C@@H](O)c1ccc(O)c(O)c1", "(S)-adrenaline"),
        ("CNCC(O)C1=CC(O)=C(O)C=C1", "4-[1-hydroxy-2-(methylamino)ethyl]benzene-1,2-diol"),
        ("NC[C@H](O)c1ccc(O)c(O)c1", "(R)-noradrenaline"),
        ("CC(C)NC[C@H](O)c1ccc(O)c(O)c1", "L-isoprenaline"),
        ("NCC(O)c1ccc(O)c(O)c1", "noradrenaline"),
        ("OC1=C(O)C=C(CCN)C=C1O", "5-Hydroxydopamine"),
        ("N[C@@H](CSc1cc(C[C@H](N)C(O)=O)cc(O)c1O)C(O)=O", "Cysteinyldopa"),
        ("CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1", "dobutamine"),
        ("NCCc1ccc(O)c(O)c1", "dopamine"),
        ("CCCN(CCC)CCC1=CC(=C(C=C1)O)O", "4-[2-(dipropylamino)ethyl]benzene-1,2-diol"),
        ("C[N+](C)(C)CCc1ccc(O)c(O)c1", "Coryneine"),
        
        # False positives (molecules that have similar fragments but not connected as a catecholamine)
        ("CCOC1=C(C=C(C=C1)CCNC(=O)CC2=CC=CC=C2)OCC", "N-[2-(3,4-diethoxyphenyl)ethyl]-2-phenylacetamide"),
        # More false positive examples from the provided data can be tested.
    ]
    
    for smi, name in test_molecules:
        result, reason = is_catecholamine(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")