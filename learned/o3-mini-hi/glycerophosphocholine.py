"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
Definition: The glycerol phosphate ester of a phosphocholine.
A true glycerophosphocholine (GPC) should contain (i) a phosphocholine headgroup – an oxygen‐bonded 
ethylene trimethylammonium fragment attached to a phosphate group – and (ii) this phosphate must be
linked to a glycerol (or glycerol‐like) moiety. (That helps to distinguish sphingomyelins or other false positives.) 

Note: This is only a heuristic method.
"""

from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a glycerophosphocholine.
    The algorithm works in two steps:
      1. It ensures that the molecule has at least one phosphorus atom and at least one of several
         possible phosphocholine headgroup SMARTS patterns. These patterns look for the fragment
         "COP(=O)(O)OCC[N+](C)(C)C" (with small variations to allow for alternative protonation states).
      2. It checks that for at least one phosphorus atom the phosphate is attached (via an oxygen)
         to a carbon whose connectivity is consistent with being derived from a glycerol (or glycerol‐like)
         backbone. In glycerophosphocholine (and phosphatidylcholines) the phosphate is attached at the
         sn‐3 (primary) position of the glycerol and that carbon tends to have few heavy‐atom neighbors.
         (This also helps exclude sphingomyelins and other phospholipids that “accidentally” contain
         a choline fragment.)
         
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      bool: True if the molecule is accepted as glycerophosphocholine, else False.
      str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1. Require at least one phosphorus atom.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found; not a phospholipid"

    # Define a few SMARTS patterns for a phosphocholine headgroup.
    # These target a phosphate bonded to a carbon which is then bonded to an ethylene fragment ending
    # in a trimethylammonium group. (They allow for variation in protonation.)
    phosphocholine_smarts_list = [
        "COP(=O)(O)OCC[N+](C)(C)C",      # all oxygens protonated
        "COP(=O)([O-])OCC[N+](C)(C)C",    # one oxygen deprotonated
        "COP([O-])(=O)OCC[N+](C)(C)C"     # alternative ordering of deprotonated oxygen
    ]
    headgroup_found = False
    for smarts in phosphocholine_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            headgroup_found = True
            break
    if not headgroup_found:
        return False, "Phosphocholine headgroup substructure not found"

    # Step 2. Look for evidence of the glycerol (or glycerol‐like) linkage.
    # We iterate over all phosphorus atoms and (for each) examine oxygen neighbors.
    # Excluding the oxygen that –by itself– carries the choline chain, we require that the phosphate
    # is connected (via an oxygen) to a carbon that is not highly substituted (heuristically, having
    # three or fewer heavy-atom bonds). This is a proxy for a glycerol (3‐carbon) backbone.
    glycerol_found = False
    for p_atom in mol.GetAtoms():
        if p_atom.GetAtomicNum() != 15:
            continue
        for nb in p_atom.GetNeighbors():
            # Look only at oxygen neighbors.
            if nb.GetAtomicNum() != 8:
                continue
            # Get oxygen neighbors excluding the phosphorus.
            oxy_neighbors = [a for a in nb.GetNeighbors() if a.GetIdx() != p_atom.GetIdx()]
            # Heuristic: if any neighbor of this oxygen is a carbon that itself is not overly substituted,
            # then assume it comes from a glycerol backbone.
            # (Also, if that carbon is connected to a positively charged nitrogen in a trimethylammonium group,
            # we assume this oxygen is already assigned to the headgroup and skip it.)
            is_choline_oxygen = False
            for onbr in oxy_neighbors:
                # Look for a nitrogen with a positive charge (i.e. part of the choline)
                if onbr.GetAtomicNum() == 7 and onbr.GetFormalCharge() == 1:
                    is_choline_oxygen = True
                    break
            if is_choline_oxygen:
                continue  # skip this oxygen because it is used for the choline branch

            # Now, check if any of the remaining neighbors is a carbon with low substitution.
            for onbr in oxy_neighbors:
                if onbr.GetAtomicNum() == 6:  # carbon
                    # Count heavy-atom (atomic num > 1) neighbors of the carbon.
                    heavy_neighbors = [a for a in onbr.GetNeighbors() if a.GetAtomicNum() > 1]
                    if len(heavy_neighbors) <= 3:
                        glycerol_found = True
                        break
            if glycerol_found:
                break
        if glycerol_found:
            break

    if not glycerol_found:
        return False, "Glycerol backbone connectivity not found (phosphate not linked to a short carbon chain)"
    
    return True, "Molecule contains a phosphocholine moiety (glycerol phosphate ester of phosphocholine) with proper glycerol connectivity"


# Example usage.
if __name__ == "__main__":
    # Test with one of the provided SMILES strings.
    test_smiles = "P(OC[C@@H](CO)OC(CCCCCCC/C=C\\CC(CCCCCC)O)=O)(=O)(OCC[N+](C)(C)C)[O-]"
    result, reason = is_glycerophosphocholine(test_smiles)
    print(result, reason)