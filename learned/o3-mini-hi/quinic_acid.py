"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: A cyclitol carboxylic acid (quinic acid and derivatives)
Quinic acid and its natural derivatives have a very characteristic cyclitol core:
  – A saturated, non‐aromatic 6-membered (cyclohexane) ring.
  – The ring carries a carboxyl group (which may be free or esterified).
  – Several oxygen substituents (e.g. hydroxy groups, acylated –OH groups) are attached.
We require that the candidate ring has at least 3 oxygen–containing substituents with 
at least one being a (esterified or free) carboxyl group.
"""

from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative based on its SMILES string.
    
    For our purposes quinic acid (and derivatives) must have:
      - A non-aromatic, saturated cyclohexane ring.
      - At least one carboxyl (or esterified carboxyl) group attached directly to the ring.
      - At least three oxygen substituents attached to the ring (including the carboxyl).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a quinic acid derivative, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information: list of rings (as tuples of atom indices)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Pre-calculate acid/ester pattern: this pattern catches the basic unit C(=O)O,
    # which is common to both free carboxylic acids and many esters.
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    # Check each candidate ring: we want a 6-membered ring that is non-aromatic, and all atoms are carbon.
    for ring in rings:
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is a non-aromatic carbon (sp3)
        candidate = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                candidate = False
                break
        if not candidate:
            continue
        
        # For this candidate cyclohexane ring, count the number of oxygen substituents
        # attached to its atoms. (Neighbors not in the ring.)
        oxy_substituents = 0
        carboxyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider neighbors not in the ring
                if nbr.GetIdx() in ring:
                    continue
                # Look for an oxygen directly attached to the ring (as -OH or as part of an ester group)
                if nbr.GetAtomicNum() == 8:
                    oxy_substituents += 1
                # Also, sometimes the carboxyl is attached via a carbon.
                # Check if a non-ring neighbor is a carbon that is part of a C(=O)O moiety.
                if nbr.GetAtomicNum() == 6:
                    # Check if this carbon is part of any match of acid_pattern.
                    # The pattern's first atom (the carbon) should match.
                    for match in acid_matches:
                        # If the carbon (match[0]) is the neighbor, then mark carboxyl found.
                        if match[0] == nbr.GetIdx():
                            carboxyl_found = True
                            break
            # In case already found carboxyl, no need to check further for that atom.
            # (We still want to count all oxygen substituents.)
        
        # We now decide on a candidate ring: must have at least one carboxyl group attached, 
        # and at least 3 oxygen substituents in total (this threshold may help suppress false positives).
        if carboxyl_found and oxy_substituents >= 3:
            return True, "Found cyclohexane ring with sufficient oxygen substituents and an acid/ester group (quinic acid core)"
    
    # No candidate ring meeting the criteria was found.
    return False, "No cyclohexane ring with required oxygen substituents and acid group detected"

# Uncomment for testing:
# test_smiles = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"  # (+)-quinic acid
# result, reason = is_quinic_acid(test_smiles)
# print(result, reason)