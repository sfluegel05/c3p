"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
#!/usr/bin/env python3
"""
Classifies: a 3-oxo-Delta(4) steroid
Definition: A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
This code checks for a tetracyclic (steroid) nucleus and for an enone motif
(i.e. a ketone conjugated with an adjacent C=C in a six-membered ring).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined as a steroid (typically with a tetracyclic nucleus)
    having a 3-oxo group (ketone) conjugated with an adjacent C=C double bond in the A ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a steroid-like nucleus by verifying that the molecule contains four rings.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4:
        return False, f"Not enough rings ({n_rings} found); expected a tetracyclic steroid nucleus"
    
    # Define the SMARTS pattern for the enone functionality in a six‐membered ring.
    # This pattern demands:
    #    - A carbon atom (in a ring of size 6) bonded to a ketonic carbon
    #    - The ketonic carbon ([CX3]) in a ring (r6) with an oxygen double bond (=O)
    #    - Followed by a carbon in the ring that is doubly-bonded to another ring carbon.
    enone_smarts = "[#6;r6][CX3;r6](=O)[#6;r6]=[#6;r6]"
    enone = Chem.MolFromSmarts(enone_smarts)
    if enone is None:
        return False, "Error in generating the enone SMARTS pattern"
    
    # Check if the molecule contains the enone substructure
    if not mol.HasSubstructMatch(enone):
        return False, "Missing the alpha,beta-unsaturated ketone motif (3-oxo & Δ(4)) in a six-membered ring"
    
    return True, "Contains a 3-oxo steroid motif with a Δ(4) unsaturation"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)O"  # (20S)-20-hydroxypregn-4-en-3-one
    result, reason = is_3_oxo_Delta_4__steroid(test_smiles)
    print(result, reason)