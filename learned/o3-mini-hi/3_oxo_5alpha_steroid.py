"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
#!/usr/bin/env python
"""
Classifies: 3-oxo-5α-steroid
Definition: A 3-oxo steroid that has alpha configuration at position 5.
Heuristic: The molecule should have a steroid-like fused ring system
(three six-membered rings and one five-membered ring), a ketone (C=O) group
that is part of one of those rings (as in a 3-oxo group) and at least one chiral
center in the steroid backbone showing “alpha” configuration (here indicated by the '@@'
stereochemistry flag).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5α-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-oxo-5α-steroid, else False.
        str: Reason explaining the result.
    """
    # Parse the SMILES string to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a steroid-like fused ring system.
    # A typical steroid has 4 rings: three six-membered rings and one five-membered ring.
    ring_info = mol.GetRingInfo().AtomRings()
    six_membered = sum(1 for ring in ring_info if len(ring) == 6)
    five_membered = sum(1 for ring in ring_info if len(ring) == 5)
    if six_membered < 3 or five_membered < 1:
        return False, "Molecule does not have the typical steroid ring system (expect three 6-membered and one 5-membered ring)"
    
    # Check for the 3-oxo (ketone) group in a ring.
    # We search for a carbonyl group on a ring carbon. The SMARTS "[#6;R]=O" matches a carbon in a ring double-bonded to oxygen.
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if len(ketone_matches) == 0:
        return False, "No ketone group (C=O in a ring) found for the 3-oxo functionality"
    
    # Further filter: we want the carbonyl carbon to be attached to at least two ring atoms.
    is_3oxo = False
    for match in ketone_matches:
        carbonyl_idx = match[0]
        atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Count if at least two neighbors are part of a ring.
        ring_neighbors = sum(1 for nb in atom.GetNeighbors() if nb.IsInRing())
        if ring_neighbors >= 2:
            is_3oxo = True
            break
    if not is_3oxo:
        return False, "Ketone group not positioned as expected in a 3-oxo steroid"
    
    # Check for the steroid backbone.
    # Here we use a simple SMARTS pattern that matches a fused polycyclic system resembling a steroid nucleus.
    # Note: This is a simplified pattern and may not capture all steroid variants.
    steroid_scaffold = Chem.MolFromSmarts("C1CC2CCC3C(C2)CCC13")
    if not mol.HasSubstructMatch(steroid_scaffold):
        return False, "Steroid backbone not detected"
    
    # Check for the alpha configuration at position 5.
    # As a heuristic, we examine atoms in the steroid scaffold: at least one chiral center specified
    # with '@@' (which is often used to indicate alpha stereochemistry as defined by CIP rules)
    scaffold_matches = mol.GetSubstructMatches(steroid_scaffold)
    alpha_found = False
    for match in scaffold_matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider atoms with defined chirality.
            if atom.HasProp('_ChiralityPossible'):
                chiral_tag = atom.GetChiralTag()
                # In RDKit, the CHI_TETRAHEDRAL_CCW tag (counterclockwise) is often encountered
                # and can be used heuristically as alpha configuration, though a more rigorous
                # assignment might be necessary in a full implementation.
                if chiral_tag == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    alpha_found = True
                    break
        if alpha_found:
            break
    if not alpha_found:
        return False, "No chiral center with the expected alpha configuration at position 5 detected"
    
    return True, "Molecule classified as a 3-oxo-5α-steroid (steroid backbone, 3-oxo group, and alpha configuration found)"

# The module can be tested with examples by calling is_3_oxo_5alpha_steroid(smiles)