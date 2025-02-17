"""
Classifies: CHEBI:22315 alkaloid
"""
#!/usr/bin/env python3
"""
Classifies: Alkaloid

Improved heuristic for classifying alkaloids:
  • Reject molecules with a very high heteroatom ratio (>0.5).
  • We count amide bonds (pattern C(=O)N) but allow a few if the molecule is complex.
  • Require at least one nitrogen atom.
  • If any nitrogen is in a ring then check that the ring system is not just a single small (≤6 atoms) fully aromatic ring.
  • If all nitrogen atoms are exocyclic, then try to identify benzylamine-like motifs.
       • We check BOTH a direct benzylamine (N–CH2–aromatic) and a relaxed pattern (N–C–CH2–aromatic).
       • Failing that, we check if any nitrogen is “close” (within 3 bonds) to a sizable (≥6 atom aromatic) ring.
  • Also, if the molecule is overall complex (many heavy atoms and rings), we lean toward an alkaloid.
"""

from rdkit import Chem
import numpy as np

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string using improved heuristics.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as an alkaloid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic properties: total atoms, heavy atoms, heteroatoms ratio.
    total_atoms = mol.GetNumAtoms()
    heavy_atom_count = mol.GetNumHeavyAtoms()
    hetero_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1,6)]
    hetero_ratio = len(hetero_atoms) / total_atoms if total_atoms > 0 else 0
    if hetero_ratio > 0.5:
        return False, "High heteroatom ratio suggests peptide/nucleotide rather than alkaloid"
    
    # Count amide bonds using pattern "C(=O)N"
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern) if amide_pattern is not None else []
    n_amide = len(amide_matches)
    # Allow multiple amide bonds if molecule is sufficiently complex.
    if n_amide > 1 and (heavy_atom_count < 25 or len(mol.GetRingInfo().AtomRings()) < 2):
        return False, "Multiple amide bonds in a small/simple molecule suggest peptide/nucleotide rather than alkaloid"
    
    # Ensure at least one nitrogen atom exists.
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found; alkaloids require at least one nitrogen."
    
    # Get ring information.
    rings = list(mol.GetRingInfo().AtomRings())
    n_rings = len(rings)
    
    # Check if any nitrogen is in a ring.
    nitro_in_ring = False
    for N in nitrogen_atoms:
        idx = N.GetIdx()
        for ring in rings:
            if idx in ring:
                nitro_in_ring = True
                break
        if nitro_in_ring:
            break

    # Helper function: checks if there is any aromatic ring with at least min_size atoms.
    def has_aromatic_ring(min_size=6):
        for ring in rings:
            if len(ring) >= min_size and all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                return True
        return False

    # Helper function: whether a given atom is within max_bonds bonds from any atom in a sizable aromatic ring.
    def atom_near_aromatic(atom_idx, max_bonds=3, min_ring_size=6):
        # Using the RDKit distance matrix.
        dmat = Chem.GetDistanceMatrix(mol)
        for ring in rings:
            if len(ring) >= min_ring_size:
                for idx in ring:
                    if dmat[atom_idx][idx] <= max_bonds:
                        return True
        return False

    # Case 1: At least one nitrogen is in a ring.
    if nitro_in_ring:
        # If only one ring, and that ring is small (<=6 atoms) and fully aromatic,
        # then it might be just a simple aromatic heterocycle.
        if n_rings == 1 and len(rings[0]) <= 6:
            if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in rings[0]):
                return False, ("Only a single small, fully aromatic ring detected; "
                               "likely a simple heterocycle (amine) rather than an alkaloid")
        return True, "Nitrogen present in a ring system with a complex or multiple-ring scaffold; matches alkaloid features"
    
    # Case 2: No nitrogen lies in a ring.
    # First, check for direct benzylamine motif (N directly attached to CH2 attached to an aromatic ring)
    benzylamine_pat_direct = Chem.MolFromSmarts("[NX3;!$(*@R)]-[CH2]-c1ccccc1")
    if benzylamine_pat_direct and mol.HasSubstructMatch(benzylamine_pat_direct):
        return True, ("Exocyclic nitrogen attached directly to a benzylic CH2 linked to an aromatic ring; "
                       "matches known alkaloid features")
    
    # Check for an extended benzylamine pattern: allow one extra intervening carbon.
    benzylamine_pat_extended = Chem.MolFromSmarts("[NX3;!$(*@R)]~[C]-[CH2]-c1ccccc1")
    if benzylamine_pat_extended and mol.HasSubstructMatch(benzylamine_pat_extended):
        return True, ("Exocyclic nitrogen with a short (one extra bond) linker to an aromatic ring; "
                       "matches known alkaloid features")
    
    # Next, if any nitrogen is near a sizable aromatic ring (within 3 bonds)
    for N in nitrogen_atoms:
        if atom_near_aromatic(N.GetIdx(), max_bonds=3, min_ring_size=6):
            # Also require a minimal molecular complexity to avoid oversimplistic amines.
            if heavy_atom_count >= 12:
                return True, ("Nitrogen found within 3 bonds of a sizable aromatic ring and molecule complexity suffices; "
                               "matches alkaloid features")
    
    # Finally, if the molecule is relatively large and complex (many heavy atoms and multiple rings)
    if heavy_atom_count >= 20 and n_rings >= 2:
        return True, ("Complex molecular scaffold with multiple rings and sufficient heavy atoms; "
                       "despite exocyclic nitrogen, matches features of alkaloids")
    
    return False, ("All nitrogen atoms are exocyclic and no benzylamine or proximity to a sizable aromatic ring was detected, "
                   "and the molecule is not complex enough; likely a simple amine rather than an alkaloid")

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
            "C[C@H](N)[C@@H](O)C1=CC=CC=C1"
    }
    
    for name, smi in test_examples.items():
        is_alka, reason = is_alkaloid(smi)
        print(f"\nMolecule: {name}\n SMILES: {smi}\n Classified as alkaloid? {is_alka}\n Reason: {reason}")