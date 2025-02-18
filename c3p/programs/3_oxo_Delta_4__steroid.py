"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-Delta(4) steroid
Definition: A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
This code checks for a fused tetracyclic (steroid) nucleus and for an enone motif 
(i.e. a ketone conjugated with an adjacent C=C in a six‐membered ring).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined as a steroid (with a fused tetracyclic nucleus) 
    having a 3-oxo group (ketone) conjugated to an adjacent C=C double bond in a six-membered ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a steroid nucleus.
    # A common steroid core is a fused tetracyclic system.
    # This SMARTS pattern is a rough representation of the steroid nucleus 
    # (three six-membered rings and one five-membered ring fused together).
    steroid_core_smarts = "C1[C@@H]2CC[C@H]3C4CC[C@H](C4)[C@H]3[C@H]2C1"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Error in generating steroid nucleus SMARTS"
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid nucleus not found"
    
    # 2. Check for the enone functionality (alpha,beta-unsaturated ketone) 
    # in a six-membered ring.
    # The SMARTS below requires:
    #   - a ring carbon ([#6;R]) bonded to a ketonic carbon ([CX3;R](=O))
    #   - that ketonic carbon bonded to another ring carbon engaged in a double bond.
    enone_smarts = "[#6;R][CX3;R](=O)[#6;R]=[#6;R]"
    enone = Chem.MolFromSmarts(enone_smarts)
    if enone is None:
        return False, "Error in generating enone SMARTS pattern"
    
    enone_matches = mol.GetSubstructMatches(enone)
    if not enone_matches:
        return False, "Missing the alpha,beta-unsaturated ketone motif (3-oxo & Δ(4))"
    
    # To further ensure that the enone motif is present in a six-membered ring (A-ring),
    # iterate over the enone matches and check if all atoms in the match are contained
    # in any one six-membered ring.
    ring_info = mol.GetRingInfo()
    enone_in_ring = False
    for match in enone_matches:
        # Iterate over every ring in the molecule.
        for ring in ring_info.AtomRings():
            if len(ring) == 6 and set(match).issubset(ring):
                enone_in_ring = True
                break
        if enone_in_ring:
            break
    if not enone_in_ring:
        return False, "Enone motif not found in a six-membered ring"
    
    return True, "Contains a 3-oxo steroid motif with a Δ(4) unsaturation"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)O"  # (20S)-20-hydroxypregn-4-en-3-one
    result, reason = is_3_oxo_Delta_4__steroid(test_smiles)
    print(result, reason)