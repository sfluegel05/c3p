"""
Classifies: CHEBI:22315 alkaloid
"""
#!/usr/bin/env python3
"""
Classifies: Alkaloid

Definition (enhanced heuristic):
  Alkaloids are naturally occurring basic nitrogen compounds that typically have at least one 
  heterocyclic ring or a benzylamine-type motif. In addition, many alkaloids feature fused or multi‐ring 
  scaffolds. In contrast, compounds whose nitrogen is exclusively exocyclic (e.g. dopamine, mescaline, 
  serotonin) are usually amines. Compounds with a very high heteroatom ratio (as seen in peptides or nucleotides) 
  are also excluded.
  
This implementation:
  1. Checks that the input SMILES is valid.
  2. Computes a heteroatom ratio to help filter out non-natural product scaffolds.
  3. Searches for nitrogen atoms in rings. However, if the only ring is a small fully aromatic heterocycle, 
     it is considered a likely false positive.
  4. Detects a benzylamine-like motif – now allowing either CH2 or CH as the benzylic linker.
  5. Uses overall molecular complexity (i.e. heavy atom count and ring count, especially fused rings) 
     as an extra cue to catch alkaloids even when nitrogen is exocyclic.
"""

from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Classifies whether a molecule is an alkaloid based on its SMILES.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        bool: True if molecule is classified as alkaloid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Compute a rough heteroatom ratio (non-carbon, non-hydrogen atoms).
    total_atoms = mol.GetNumAtoms()
    hetero_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)]
    hetero_ratio = len(hetero_atoms) / total_atoms if total_atoms > 0 else 0
    if hetero_ratio > 0.5:
        return False, "High heteroatom ratio suggests peptide/nucleotide rather than alkaloid"
    
    # Ensure there is at least one nitrogen.
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found; alkaloids require at least one nitrogen."
    
    # Retrieve ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    ring_list = list(ring_info)
    n_rings = len(ring_list)
    
    # Check for fused rings (i.e. rings that share at least 2 atoms)
    fused = False
    for i in range(len(ring_list)):
        for j in range(i+1, len(ring_list)):
            if len(set(ring_list[i]).intersection(ring_list[j])) >= 2:
                fused = True
                break
        if fused:
            break

    # Determine if any nitrogen is part of a ring.
    nitrogen_in_ring = False
    for atom in nitrogen_atoms:
        idx = atom.GetIdx()
        for ring in ring_list:
            if idx in ring:
                nitrogen_in_ring = True
                break
        if nitrogen_in_ring:
            break

    # If a nitrogen lies in a ring, make further checks.
    if nitrogen_in_ring:
        # If there is exactly one ring and it is small (6 atoms or fewer),
        # check if it is fully aromatic. In that case, it is likely just a simple heterocycle.
        if n_rings == 1 and len(ring_list[0]) <= 6:
            if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring_list[0]):
                return False, ("Only a single small, fully aromatic ring detected; "
                               "likely a simple heterocycle (amine) rather than an alkaloid")
        # Otherwise, if the ring system is complex (multiple or fused rings), we take it as a match.
        return True, "Nitrogen present in a ring system with a complex scaffold; matches features of alkaloids"
    
    # Case where no nitrogen is directly in a ring:
    # Try to detect a benzylamine-like motif.
    # Pattern 1: exocyclic nitrogen attached to a CH2 group that is connected to an aromatic ring.
    benzylamine_pattern1 = Chem.MolFromSmarts("[NX3;!$(*@R)]-[CH2]-c1ccccc1")
    # Pattern 2: exocyclic nitrogen attached to a benzylic CH (which may be chiral) linked to an aromatic ring.
    benzylamine_pattern2 = Chem.MolFromSmarts("[NX3;!$(*@R)]-[CH]-c1ccccc1")
    if (benzylamine_pattern1 is not None and mol.HasSubstructMatch(benzylamine_pattern1)) or \
       (benzylamine_pattern2 is not None and mol.HasSubstructMatch(benzylamine_pattern2)):
        return True, ("Exocyclic nitrogen attached to a benzylic (CH2 or CH) center linked to an aromatic ring; "
                       "matches known alkaloid features")
    
    # Finally, if the molecule has a display of overall scaffold complexity (many heavy atoms and multiple rings)
    # then classify it as an alkaloid even if no nitrogen is inside a ring.
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count >= 20 and n_rings >= 2:
        return True, ("Complex molecular scaffold with multiple rings and sufficient heavy atoms; "
                       "despite exocyclic nitrogen, matches features of alkaloids")
    
    return False, ("All nitrogen atoms are exocyclic, no benzylamine motif detected, and the scaffold is not complex; "
                   "likely an amine rather than an alkaloid")


# Example test cases to illustrate the classification outcomes:
if __name__ == "__main__":
    test_examples = {
        "(-)-selegiline": "C[C@H](Cc1ccccc1)N(C)CC#C",
        "Protoveratrine B": "CC[C@@H](C)C(=O)O[C@H]1[C@H](O)[C@H]2[C@@H](CN3C[C@@H](C)CC[C@H]3[C@@]2(C)O)[C@@H]2C[C@@]34O[C@@]5(O)[C@@H]([C@@H](OC(C)=O)[C@@H](OC(C)=O)[C@H]3[C@]12O)[C@]4(C)CC[C@@H]5OC(=O)[C@@](C)(O)[C@@H](C)O",
        "Thalmine": "COc1ccc2C[C@@H]3N(C)CCc4cc(OC)c(Oc5c(O)c(OC)cc6[C@H](Cc7ccc(Oc1c2)cc7)N(C)CCc56)cc34"
    }
    for name, smi in test_examples.items():
        result, reason = is_alkaloid(smi)
        print(f"\nMolecule: {name}\n SMILES: {smi}\n Classified as alkaloid? {result}\n Reason: {reason}")