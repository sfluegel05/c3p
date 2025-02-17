"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar – any sugar having one or more alcoholic hydroxy groups replaced by 
substituted or unsubstituted amino groups.

The revised algorithm:
  1. Parse the molecule from the SMILES.
  2. Loop over all rings (from RDKit) that are 5 or 6 members.
  3. In each such ring, require that (a) it contains at least one oxygen (typical for a sugar)
     and (b) most other atoms are carbons.
  4. For each ring atom that is carbon, examine its neighbors not in the ring.
       • For oxygen atoms, count them if they are –OH (i.e. have at least one hydrogen).
       • For nitrogen atoms, only count them if they are sp3 (not aromatic) and if their local 
         environment suggests they are amino substituents. In particular we mark as “acceptable”
         a nitrogen if either it carries a free –NH2 (at least one attached hydrogen) 
         or it is part of an acetamido group (we try to match a carbonyl–methyl via SMARTS).
  5. Only if the ring overall has at least three exocyclic substituents and at least one acceptable
     amino substituent do we classify the molecule as an amino sugar.
     
If a ring is found but no acceptable amino substituent meets the criteria then we return False 
with a reason indicating that sugar-like rings were detected but no “amino” group was found.
If no sugar-like ring is found at all, we also report that.

Because this heuristic is not perfect the program may (as before) mis‐classify some edge cases.
"""

from rdkit import Chem

# Pre-compile a SMARTS for an acetamido environment.
# This SMARTS looks for a nitrogen bound to a carbonyl carbon with a methyl group.
acetamido_smarts = Chem.MolFromSmarts("[NX3;!R][C](=O)[CH3]")
    
def is_valid_amino(n_atom):
    """
    Checks whether the given nitrogen atom appears to be a valid amino substituent:
       - It is sp3 hybridized (i.e. not aromatic)
       - It either has one or more hydrogen atoms (free amine) 
         or is part of an acetamido group (the nitrogen is attached to a C(=O)CH3).
    """
    # Only consider sp3 non-aromatic nitrogen
    if n_atom.GetAtomicNum() != 7:
        return False
    if n_atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False
    # Check if it has at least one attached hydrogen.
    if n_atom.GetTotalNumHs() > 0:
        return True
    # Otherwise see if its neighborhood matches the acetamido pattern.
    # Create a temporary molecule consisting only of the neighbor environment.
    env = Chem.PathToSubmol(n_atom.GetOwningMol(), [b.GetIdx() for b in n_atom.GetBonds()])
    # If the acetamido pattern is found anywhere return True.
    if env.HasSubstructMatch(acetamido_smarts):
        return True
    return False

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    
    The algorithm inspects rings of size 5 or 6. For each candidate ring:
       • The ring must contain at least one oxygen (and ideally the remaining atoms are carbons).
       • Exocyclic substituents attached to ring carbons are examined.
             - Hydroxy groups (oxygen with hydrogen) count as substituents.
             - Nitrogen substituents are only accepted if they are sp3 (non‐aromatic)
               and if they appear to be free amines or part of an acetamido moiety.
       • The candidate ring must have at least three exocyclic substituents overall,
         and at least one acceptable amino substituent.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       (bool, str): True (with a reason message) if the molecule is classified as an amino sugar;
                    otherwise, False with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    candidate_found = False

    for ring in ring_info:
        # Consider only rings of length 5 or 6.
        if len(ring) not in (5, 6):
            continue
        
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count how many oxygen atoms are in the ring.
        oxy_in_ring = [a for a in ring_atoms if a.GetAtomicNum() == 8]
        # In typical sugar rings we expect at least one oxygen; sometimes more may occur,
        # but we demand that the majority of ring atoms should be either carbon or oxygen.
        if len(oxy_in_ring) < 1:
            continue  # not sugar-like
        
        # Optionally, we may prefer rings that have one oxygen (pyranose or furanose rings)
        # but here we are slightly more forgiving.
        
        # Examine exocyclic substituents on ring carbons.
        total_exo = 0
        amino_found = False
        for atom in ring_atoms:
            # We mainly check exocyclic substituents on carbons.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms in the ring
                at_num = nbr.GetAtomicNum()
                # For an oxygen substituent, count if it likely is an -OH (check for attached Hs)
                if at_num == 8:
                    # Count only if there is at least one hydrogen (note: explicit or implicit)
                    if nbr.GetTotalNumHs() > 0:
                        total_exo += 1
                    else:
                        total_exo += 1  # even if not explicit, count it as a substituent
                elif at_num == 7:
                    # For nitrogen, check if it is a valid amino substituent.
                    if is_valid_amino(nbr):
                        total_exo += 1
                        amino_found = True
        # Require at least 3 exocyclic substituents to be “sugar‐like.”
        if total_exo < 3:
            continue
        candidate_found = True
        if amino_found:
            return True, ("Sugar‐like ring (5- or 6-membered, with at least one oxygen and at least three exocyclic substituents) "
                          "detected that has an acceptable amino substituent (free –NH2 or acetamido group).")
    # End ring loop.
    if not candidate_found:
        return False, "No typical sugar ring (5- or 6-membered ring with at least one oxygen and multiple exocyclic substituents) was detected."
    else:
        return False, ("Sugar-like ring(s) were detected, but none contained an acceptable amino substituent "
                       "(a free amino or acetamido group); hence the molecule was not classified as an amino sugar.")

# Example usage:
if __name__ == "__main__":
    # Example SMILES that should be classified as an amino sugar by our improved algorithm.
    test_smiles = ("O([C@H]1[C@@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2NC(=O)C)CO)"
                   "[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO")
    result, reason = is_amino_sugar(test_smiles)
    print("Result:", result)
    print("Reason:", reason)