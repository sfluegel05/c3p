"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
#!/usr/bin/env python3
"""
Classifies: Glycosaminoglycan (defined as any polysaccharide containing 
a substantial proportion of aminomonosaccharide residues)
"""

from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    We use the heuristic that a glycosaminoglycan is a polysaccharide (has at least 
    two sugar rings) in which at least half of the candidate sugar rings are 
    aminomonosaccharide units (i.e. one of the ring carbons carries a substituent 
    nitrogen).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a glycosaminoglycan, False otherwise.
        str: Reason for the classification decision.
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information using RDKit's ring info
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # list of tuples with atom indices for each ring
    
    # Heuristic: Look for candidate sugar rings
    # For a pyranose, we expect a 6-membered ring with exactly one oxygen (usually the ring heteroatom).
    candidate_sugar_rings = []
    aminoglycose_flags = []  # parallel list: True if ring has N substituent on one of its carbons
    
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # only consider 6-membered rings
        # Count atoms in ring that are oxygen. (Our candidate sugar ring should have exactly one O in the ring.)
        oxygen_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxygen_in_ring != 1:
            continue
        
        # Consider this ring a candidate sugar unit.
        candidate_sugar_rings.append(ring)
        
        # Now check if any carbon in the ring is substituted with a nitrogen outside the ring.
        amino_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We only check carbons (atomic number 6) in the ring
            if atom.GetAtomicNum() != 6:
                continue
            for neighbor in atom.GetNeighbors():
                # Only consider neighbors that are not in the ring
                if neighbor.GetIdx() in ring:
                    continue
                # Check if the substituent is a nitrogen (atomic number 7)
                if neighbor.GetAtomicNum() == 7:
                    amino_found = True
                    break
            if amino_found:
                break
        aminoglycose_flags.append(amino_found)
    
    total_sugars = len(candidate_sugar_rings)
    if total_sugars < 2:
        return False, f"Not a polysaccharide: found only {total_sugars} candidate sugar ring(s)."
    
    amino_count = sum(1 for flag in aminoglycose_flags if flag)
    ratio = amino_count / total_sugars
    
    if ratio >= 0.5:
        return True, f"Polysaccharide with {total_sugars} sugar rings, of which {amino_count} ({ratio*100:.1f}%) are amino sugars."
    else:
        return False, f"Candidate polysaccharide with {total_sugars} sugar rings but only {amino_count} ({ratio*100:.1f}%) are amino sugars."
        
# Example usage:
if __name__ == "__main__":
    # Test with one of the provided SMILES (this is just one example, not necessarily a perfect glycosaminoglycan)
    test_smiles = "O=C1O[C@@H](C=C[C@@H](CC[C@]23O[C@H](C=4C(=C(NC(CC=C1C)=O)C=C(O)C4)O2)[C@@H](C)C(C3)=O)CC)[C@@H](O)C=C(C)C"
    result, reason = is_glycosaminoglycan(test_smiles)
    print("Result:", result)
    print("Reason:", reason)