"""
Classifies: CHEBI:22315 alkaloid
"""
#!/usr/bin/env python3
"""
Classifies: Alkaloid

Improved heuristic:
  • Reject molecules with high heteroatom fraction (>0.5) or multiple amide bonds (likely peptides)
  • Require at least one nitrogen atom.
  • If any nitrogen is in a ring then check that the ring system is not just a single small fully aromatic ring.
  • If no nitrogen is inside a ring, then either detect a benzylamine motif (allowing CH2 or CH as benzylic linker)
    or check that at least one nitrogen is “close” (within two bonds) to a sizable aromatic ring.
  • Also require a minimum molecular complexity (by heavy atom count) when only exocyclic nitrogens are present.
"""

from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determine whether a molecule is an alkaloid based on its SMILES string using improved heuristics.
    
    Args:
        smiles (str): Input SMILES string.
    
    Returns:
        bool: True if molecule is classified as an alkaloid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate heteroatom ratio (atoms other than C and H)
    total_atoms = mol.GetNumAtoms()
    hetero_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() not in (1,6)]
    hetero_ratio = len(hetero_atoms) / total_atoms if total_atoms > 0 else 0
    if hetero_ratio > 0.5:
        return False, "High heteroatom ratio suggests peptide/nucleotide rather than alkaloid"
    
    # Count amide bonds (pattern C(=O)N) to weed out peptides.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern) if amide_pattern is not None else []
    if len(amide_matches) > 1:
        return False, "Multiple amide bonds suggest a peptide/nucleotide rather than a typical alkaloid"
    
    # Ensure at least one nitrogen is present.
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found; alkaloids require at least one nitrogen."
    
    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    rings = list(ring_info)
    n_rings = len(rings)
    
    # Check if any nitrogen atom is in a ring.
    nitro_in_ring = False
    for N in nitrogen_atoms:
        idx = N.GetIdx()
        for ring in rings:
            if idx in ring:
                nitro_in_ring = True
                break
        if nitro_in_ring:
            break

    # Helper: check if there is any aromatic ring of at least given size.
    def has_aromatic_ring(min_size=6):
        for ring in rings:
            if len(ring) >= min_size and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                return True
        return False

    # Helper: check if a given atom index is within 'max_bonds' of any atom in a sizable aromatic ring.
    def atom_near_aromatic(atom_idx, max_bonds=2, min_ring_size=6):
        # Get distances from atom_idx to all other atoms.
        dists = Chem.GetDistanceMatrix(mol)
        for ring in rings:
            if len(ring) >= min_ring_size:
                for idx in ring:
                    if dists[atom_idx][idx] <= max_bonds:
                        return True
        return False

    heavy_atom_count = mol.GetNumHeavyAtoms()
    
    # Case 1: A nitrogen is directly in a ring.
    if nitro_in_ring:
        # If only one ring and it is small (<=6 atoms) and fully aromatic,
        # then it is most likely just a simple heterocycle.
        if n_rings == 1 and len(rings[0]) <= 6:
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in rings[0]):
                return False, ("Only a single small, fully aromatic ring detected; "
                               "likely a simple heterocycle (amine) rather than an alkaloid")
        return True, "Nitrogen present in a ring system with a complex or multiple-ring scaffold; matches features of alkaloids"
    
    # Case 2: No nitrogen is in a ring.
    # Check for benzylamine type motifs.
    benzylamine_pat1 = Chem.MolFromSmarts("[NX3;!$(*@R)]-[CH2]-c1ccccc1")
    benzylamine_pat2 = Chem.MolFromSmarts("[NX3;!$(*@R)]-[CH]-c1ccccc1")
    if (benzylamine_pat1 and mol.HasSubstructMatch(benzylamine_pat1)) or \
       (benzylamine_pat2 and mol.HasSubstructMatch(benzylamine_pat2)):
        return True, ("Exocyclic nitrogen attached to a benzylic (CH2 or CH) center linked to an aromatic ring; "
                       "matches known alkaloid features")
    
    # If no direct benzylamine is found, look for a “neighborhood” criterion:
    # That is, if any nitrogen is within 2 bonds of a sizable aromatic ring.
    for N in nitrogen_atoms:
        if atom_near_aromatic(N.GetIdx(), max_bonds=2, min_ring_size=6):
            # Also require a minimal complexity (here we choose heavy_atom_count>=12).
            if heavy_atom_count >= 12:
                return True, ("Nitrogen found close (within 2 bonds) to a sizable aromatic ring "
                               "and molecule is not overly simple; matches alkaloid features")
    
    # Also, if the molecule is relatively complex (multiple rings and heavy atoms) even with exocyclic N,
    # it may be a natural-product alkaloid.
    if heavy_atom_count >= 20 and n_rings >= 2:
        return True, ("Complex molecular scaffold with multiple rings and sufficient heavy atoms; "
                       "despite exocyclic nitrogen, matches features of alkaloids")
    
    return False, ("All nitrogen atoms are exocyclic and no benzylamine or proximity to aromatic ring detected, "
                   "and molecule is not complex enough; likely a simple amine rather than an alkaloid")


# Example test cases to illustrate classification outcomes.
if __name__ == "__main__":
    test_examples = {
        "1-[(1R)-1-(hydroxymethyl)-7-methoxy-2,9-dimethyl-1'-spiro[1,3-dihydropyrido[3,4-b]indole-4,4'-piperidine]yl]-2-(4-morpholinyl)ethanone":
            "CN1CC2(CCN(CC2)C(=O)CN3CCOCC3)C4=C([C@@H]1CO)N(C5=C4C=CC(=C5)OC)C",
        "Protoveratrine B":
            "CC[C@@H](C)C(=O)O[C@H]1[C@H](O)[C@H]2[C@@H](CN3C[C@@H](C)CC[C@H]3[C@@]2(C)O)[C@@H]2C[C@@]34O[C@@]5(O)[C@@H]([C@@H](OC(C)=O)[C@@H](OC(C)=O)[C@H]3[C@]12O)[C@]4(C)CC[C@@H]5OC(=O)[C@@](C)(O)[C@@H](C)O",
        "(-)-selegiline":
            "C[C@H](Cc1ccccc1)N(C)CC#C",
        "hordenine":
            "CNCC(O)c1ccc(O)cc1",
        "norbelladine":
            "C1(=C(C=C(C=C1)CNCCC2=CC=C(C=C2)O)O)O",
        "cathine":
            "C[C@H](N)[C@@H](O)C1=CC=CC=C1",
    }
    
    for name, smi in test_examples.items():
        is_alka, reason = is_alkaloid(smi)
        print(f"\nMolecule: {name}\n SMILES: {smi}\n Classified as alkaloid? {is_alka}\n Reason: {reason}")