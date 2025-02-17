"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
"""
Classifies: Mucopolysaccharide (glycosaminoglycan)
Definition: “Any of the group of polysaccharides composed of alternating units from uronic acids
and glycosamines, and commonly partially esterified with sulfuric acid.”
Heuristic:
  1. Parse the molecule.
  2. Iterate over all rings (using RDKit’s ring info) and consider only rings of size 5 or 6 that contain exactly one oxygen atom (a common feature of sugar rings).
  3. For each such ring, examine atoms (mainly carbons) that are part of the ring. If a carbon atom carries a double‐bonded oxygen (as a substituent outside the ring), count that ring as “uronic acid–like.” If a carbon atom carries a nitrogen substituent, count that ring as “glycosamine–like.”
  4. Also search for sulfate ester groups ("OS(=O)(=O)O").
  5. If there are at least two “uronic acid–like” rings and two “glycosamine–like” rings (with counts roughly balanced, i.e. their difference is ≤1), classify the molecule as a mucopolysaccharide.
Note: This is a rough heuristic and may miss molecules or cause false negatives.
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string using a relaxed heuristic. 
    It looks for sugar-like rings (5- or 6-membered rings with one oxygen) and then checks if 
    substituents on ring carbons indicate either a carboxylic acid group (as a sign of a uronic acid unit)
    or an amino group (as a sign of a glycosamine unit). Additionally, it searches for sulfate ester groups.
    
    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      bool: True if the molecule is classified as a mucopolysaccharide, False otherwise.
      str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Setup counters.
    uronic_count = 0
    glyco_count = 0

    # Get ring information.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Iterate over rings (we focus on 5- or 6-membered rings, typical for sugars).
    for ring in ring_info:
        if len(ring) not in (5, 6):
            continue
        # Count the number of oxygen atoms within the ring.
        oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        # A typical sugar ring has exactly one ring oxygen.
        if oxy_in_ring != 1:
            continue

        # Initialize flags for this ring.
        is_uronic = False
        is_glyco = False

        # For each atom in the ring, if it is carbon we look at its neighbors outside the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # we only consider carbon atoms here as the backbone of a sugar ring
            if atom.GetAtomicNum() != 6:
                continue
            # Check neighbors not part of the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # Check for uronic acid-like: a carbon substituent with a double-bonded oxygen.
                # Look if the neighbor is oxygen AND the bond is a double bond.
                if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() >= 2:
                    is_uronic = True
                # Check for glycosamine-like: a substituent that is a nitrogen.
                if nbr.GetAtomicNum() == 7:
                    is_glyco = True
                    
        # Count the ring according to its substituents.
        if is_uronic:
            uronic_count += 1
        if is_glyco:
            glyco_count += 1

    # Additionally, match for sulfate ester groups using a simple SMARTS.
    sulfate_smarts = "OS(=O)(=O)O"
    sulfate_query = Chem.MolFromSmarts(sulfate_smarts)
    sulfate_matches = mol.GetSubstructMatches(sulfate_query) if sulfate_query is not None else []
    n_sulfate = len(sulfate_matches)

    # Build a summary message.
    summary = (f"Detected {uronic_count} uronic acid-like ring(s) and "
               f"{glyco_count} glycosamine-like ring(s).")
    if n_sulfate > 0:
        summary += f" Also found {n_sulfate} sulfate ester group(s)."
    
    # For mucopolysaccharides we require at least 2 of each type and roughly alternating (difference no more than 1).
    if uronic_count < 2 or glyco_count < 2:
        reason = ("Not enough sugar-like units: found "
                  f"{uronic_count} uronic acid-like and {glyco_count} glycosamine-like rings "
                  "(need at least 2 of each)")
        return False, reason
    
    if abs(uronic_count - glyco_count) > 1:
        reason = ("Repeating unit counts are not balanced (difference > 1): "
                  f"{uronic_count} uronic acid-like vs {glyco_count} glycosamine-like rings")
        return False, reason
    
    return True, summary