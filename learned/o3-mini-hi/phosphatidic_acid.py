"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: Phosphatidic acid (PA)
Definition:
  Phosphatidic acid is defined as a glycerol derivative in which one hydroxy group is esterified 
  with phosphoric acid and the other two are esterified with fatty acids.

Improvements in this version:
  - Identify a phosphate group by looking for a P atom with >=3 O neighbors and exactly one oxygen 
    (among these) that connects to a carbon (the phosphate ester).
  - Identify fatty acid ester groups using the substructure "OC(=O)" but filter out any cases where 
    the bridging oxygen is directly attached to a phosphorus.
  - When picking the glycerol carbon from each fatty acid ester, avoid the carbonyl carbon (which 
    should have a double bond to O) and choose the candidate that appears to be part of the glycerol.
  - Verify that the two fatty acid ester glycerol carbons plus the phosphate ester glycerol carbon form 
    a contiguous set of three carbons (i.e. pairwise graph distances are ≤2).
  - Check overall molecular weight for sanity.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    Phosphatidic acid is defined as a glycerol derivative in which one hydroxy group is esterified 
    with phosphoric acid and the other two are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a phosphatidic acid.
        str: Explanation of the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------------------------
    # Criterion 1: Identify a suitable phosphate group.
    # Look for a phosphorus (P) atom that has at least three oxygen neighbors.
    # Among these, exactly one oxygen should be bound to a carbon (the glycerol linkage).
    phosphate_found = False
    valid_P_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "P":
            continue
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"]
        if len(oxy_neighbors) < 3:
            continue
        # Count oxygens that are attached to some carbon (other than P)
        o_attached_to_C = 0
        for o in oxy_neighbors:
            for nbr2 in o.GetNeighbors():
                if nbr2.GetIdx() == atom.GetIdx():
                    continue
                if nbr2.GetSymbol() == "C":
                    o_attached_to_C += 1
                    break
        if o_attached_to_C == 1:
            phosphate_found = True
            valid_P_atom = atom
            break
    if not phosphate_found:
        return False, "No phosphate group found with >=3 oxygen neighbors and exactly one oxygen attached to a carbon"
    
    # -------------------------
    # Criterion 2: Identify fatty acid ester groups.
    # We search for substructure "OC(=O)". Then, we filter out those oxygen atoms that are
    # directly attached to any phosphorus (they would be part of the phosphate head).
    fa_pattern = Chem.MolFromSmarts("OC(=O)")
    fa_matches = mol.GetSubstructMatches(fa_pattern)
    
    fa_candidates = []  # Will collect candidate glycerol carbon indices from FA ester bonds.
    bridging_O_indices = set()  # For later reference to avoid using the same O in the phosphate part.
    
    for match in fa_matches:
        # match is a tuple of indices corresponding to O, C, O (where second atom is the carbonyl carbon).
        bridging_O_idx = match[0]
        bridging_O = mol.GetAtomWithIdx(bridging_O_idx)
        # Exclude if the bridging oxygen is directly attached to any phosphorus.
        if any(neighbor.GetSymbol() == "P" for neighbor in bridging_O.GetNeighbors()):
            continue
        # Now, from bridging O, choose a neighboring carbon that is likely attached as glycerol.
        candidate = None
        for nbr in bridging_O.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            # Distinguish the glycerol-attached C from the carbonyl carbon.
            # If the carbon has a double-bonded O, consider it as the carbonyl and skip.
            has_double_Bound_O = False
            for nbr2 in nbr.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                if nbr2.GetSymbol() == "O" and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_double_Bound_O = True
                    break
            if not has_double_Bound_O:
                candidate = nbr.GetIdx()
                break
        if candidate is not None:
            fa_candidates.append(candidate)
            bridging_O_indices.add(bridging_O_idx)
    
    if len(fa_candidates) != 2:
        return False, f"Found {len(fa_candidates)} fatty acid ester group(s) after filtering; expected exactly 2."
    
    # -------------------------
    # Criterion 3: Identify the glycerol backbone candidate from the phosphate.
    # From the valid P atom, pick the oxygen (other than the ones already used as FA bridging O)
    # that leads to a carbon. This should be the glycerol carbon attached via the phosphate ester.
    phospho_candidate_carbons = set()
    for o in valid_P_atom.GetNeighbors():
        if o.GetSymbol() != "O":
            continue
        if o.GetIdx() in bridging_O_indices:
            continue
        for nbr in o.GetNeighbors():
            if nbr.GetIdx() == valid_P_atom.GetIdx():
                continue
            if nbr.GetSymbol() == "C":
                phospho_candidate_carbons.add(nbr.GetIdx())
                break  # Assume first qualifying carbon is sufficient.
    
    if len(phospho_candidate_carbons) != 1:
        return False, "Did not find exactly one candidate glycerol carbon from the phosphate (found %d)" % len(phospho_candidate_carbons)
    
    phospho_carbon = list(phospho_candidate_carbons)[0]
    
    # Combine the candidate glycerol carbons: the two from FA esters and the one from the phosphate.
    glycerol_set = set(fa_candidates)
    glycerol_set.add(phospho_carbon)
    if len(glycerol_set) != 3:
        return False, "Fatty acid ester groups and the phosphate ester do not attach to three distinct glycerol carbons"
    
    # Verify that these three carbons are adjacent.
    # Check that the graph distance between any two is ≤2.
    dmat = Chem.GetDistanceMatrix(mol)
    for i in glycerol_set:
        for j in glycerol_set:
            if i == j:
                continue
            if dmat[i][j] > 2:
                return False, "The three candidate glycerol carbons are not sufficiently connected (graph distance >2)"
    
    # -------------------------
    # Criterion 4: Molecular weight check (sanity check).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a phosphatidic acid"
    
    return True, "Contains a phosphate group and exactly two fatty acid esters attached to a connected three-carbon glycerol backbone"

# Example usage:
if __name__ == "__main__":
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Is phosphatidic acid?", result)
    print("Reason:", reason)