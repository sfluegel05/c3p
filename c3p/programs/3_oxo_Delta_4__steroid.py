"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-Delta(4) steroid
Definition: A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
This code verifies that the molecule contains a fused tetracyclic (steroid) nucleus and that an 
alpha,beta-unsaturated ketone (enone) is present in one of its six-membered rings.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined as a steroid (with a fused tetracyclic nucleus with 
    three six-membered rings and one five-membered ring, typically totaling 17 carbons) having 
    a 3-oxo group that is conjugated to an adjacent C=C double bond (enone) in a six-membered ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------------------------
    # 1. Check for steroid nucleus
    # ---------------------------
    # We use a heuristic based on the typical steroid nucleus:
    #   - It is a fused tetracyclic system.
    #   - It must contain at least three 6-membered rings and one 5-membered ring.
    #   - Their union should cover roughly 17 carbon atoms.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    count6 = sum(1 for r in rings if len(r) == 6)
    count5 = sum(1 for r in rings if len(r) == 5)
    
    if count6 < 3 or count5 < 1:
        return False, "Steroid nucleus not found (insufficient six-membered and/or five-membered rings)"
    
    # Build a union of all atoms belonging to rings of size 5 or 6.
    steroid_atoms = set()
    for r in rings:
        if len(r) in [5,6]:
            steroid_atoms.update(r)
    # A steroid nucleus usually comprises 17 carbons;
    # here we require at least 15 heavy atoms in the fused ring system (to allow for some substitutions).
    if len(steroid_atoms) < 15:
        return False, "Steroid nucleus not found (fused ring atom count too low)"
    
    # ---------------------------
    # 2. Check for enone functionality within a six-membered ring
    # ---------------------------
    # We define an enone as a ketone (C(=O)) directly attached to a C=C double bond.
    enone_smarts = "C(=O)C=C"
    enone = Chem.MolFromSmarts(enone_smarts)
    if enone is None:
        return False, "Error in generating enone SMARTS pattern"
    
    enone_matches = mol.GetSubstructMatches(enone)
    if not enone_matches:
        return False, "Missing 3-oxo / Δ(4) enone motif"
    
    # Ensure that the enone match is part of a six-membered ring.
    enone_in_six = False
    for match in enone_matches:
        for ring in rings:
            if len(ring) == 6 and set(match).issubset(ring):
                enone_in_six = True
                break
        if enone_in_six:
            break
    if not enone_in_six:
        return False, "Enone motif not found in a six-membered ring"
    
    return True, "Contains a fused tetracyclic (steroid) nucleus with an alpha,beta-unsaturated ketone (3-oxo Δ(4))"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test example: (20S)-20-hydroxypregn-4-en-3-one
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)O"
    result, reason = is_3_oxo_Delta_4__steroid(test_smiles)
    print(result, reason)