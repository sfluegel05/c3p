"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: Phosphatidic acid (PA)
Definition:
  Phosphatidic acid is defined as a glycerol derivative in which one hydroxy group is esterified 
  with phosphoric acid and the other two are esterified with fatty acids.
Improvements in this version:
  - Identify a phosphate group with ≥3 oxygen neighbors and exactly one oxygen connected to a carbon.
  - Identify exactly two fatty acid ester groups (matching OC(=O)) while filtering out those directly linked to P.
  - Collect the glycerol carbons: the two from the fatty acid ester bonds and the one from the phosphate ester.
  - Check that these three carbons form a continuous unit (the glycerol backbone) by confirming that the graph‐distance 
    between any two of them is at most 2.
  - Ensure an overall minimum molecular weight (>300 Da).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    Phosphatidic acid is defined here as a glycerol derivative in which one hydroxy group is esterified 
    with phosphoric acid and the other two are esterified with fatty acids.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as phosphatidic acid.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------------------------
    # Criterion 1: Identify a suitable phosphate group.
    # The phosphate should be a phosphorus (P) atom that has at least three oxygen neighbors.
    # In addition, exactly one of its oxygen neighbors (the one forming the glycerol attachment)
    # should be bound to a carbon (i.e. not coming from a fatty acid or other headgroup).
    phosphate_found = False
    phospho_carbon_candidate = None  # This will be the carbon bonded to an oxygen of P (the phosphorylated –OH).
    valid_P_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "P":
            continue
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"]
        if len(oxy_neighbors) < 3:
            continue
        # Count how many oxygens are attached to a carbon (and not the P itself)
        o_attached_to_C = 0
        for o in oxy_neighbors:
            for nbr in o.GetNeighbors():
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetSymbol() == "C":
                    o_attached_to_C += 1
                    break
        if o_attached_to_C == 1:
            phosphate_found = True
            valid_P_atom = atom
            break  # take the first qualifying phosphate group
    if not phosphate_found:
        return False, "No phosphate group found with >=3 oxygen neighbors and exactly one oxygen attached to a carbon"
    
    # -------------------------
    # Criterion 2: Identify fatty acid ester groups.
    # We search for the substructure matching "OC(=O)".
    fa_pattern = Chem.MolFromSmarts("OC(=O)")
    fa_matches = mol.GetSubstructMatches(fa_pattern)
    
    # We will filter out matches in which the bridging O (first atom in match) is directly attached to any phosphorus.
    # Also, we record the carbon attached to that bridging oxygen (our candidate glycerol carbons).
    filtered_fa_matches = []
    fa_glycerol_carbons = set()  # indices of carbons that link to fatty acid esters (should be 2 in PA)
    bridging_O_indices = set()   # indices of oxygen atoms used as bridges in these matches.
    
    for match in fa_matches:
        bridging_O_idx = match[0]  # index of the oxygen in "O C(=O)"
        bridging_O = mol.GetAtomWithIdx(bridging_O_idx)
        # Exclude if the bridging oxygen is directly attached to any phosphorus (would be part of the phosphate head)
        if any(neighbor.GetSymbol() == "P" for neighbor in bridging_O.GetNeighbors()):
            continue
        filtered_fa_matches.append(match)
        bridging_O_indices.add(bridging_O_idx)
        # Now, among the neighbors of bridging_O, record the first carbon (note: there might be more than one,
        # but we take the one that is not the carbonyl carbon since that is part of the fatty acid chain).
        for nbr in bridging_O.GetNeighbors():
            if nbr.GetSymbol() == "C":
                # We try to avoid the carbon in the ester carbonyl (which typically has a double bond to O)
                # by checking if the neighbor has a double-bonded O (a carbonyl).
                has_dblO = any((nbr2.GetSymbol() == "O" and mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx()).GetBondTypeAsDouble() == 2) 
                               for nbr2 in nbr.GetNeighbors())
                if not has_dblO:
                    fa_glycerol_carbons.add(nbr.GetIdx())
                    break

    if len(filtered_fa_matches) != 2:
        return False, f"Found {len(filtered_fa_matches)} fatty acid ester group(s) after filtering; expected exactly 2."

    # -------------------------
    # Criterion 3: Identify the glycerol backbone.
    # In phosphatidic acid the two fatty acid chains (via the bridging oxygens) and the phosphate head
    # should all be attached to a glycerol backbone (3 connected carbons).
    # So far we have candidate glycerol carbons from the fatty acid esters; now we choose the one from the phosphate.
    phospho_candidate_carbons = set()
    for o in valid_P_atom.GetNeighbors():
        if o.GetSymbol() != "O":
            continue
        # Skip if this oxygen is one of the bridging oxygens already used in fatty acid groups.
        if o.GetIdx() in bridging_O_indices:
            continue
        # Check if this oxygen is bonded to a carbon and record that carbon.
        for nbr in o.GetNeighbors():
            if nbr.GetSymbol() == "C":
                phospho_candidate_carbons.add(nbr.GetIdx())
                break  # consider the first such carbon

    if len(phospho_candidate_carbons) != 1:
        return False, "Did not find exactly one candidate glycerol carbon from the phosphate (found %d)" % len(phospho_candidate_carbons)
    phospho_carbon = list(phospho_candidate_carbons)[0]

    # Now, the candidate glycerol backbone is the union of the two fatty acid linked carbons plus the phosphate-linked carbon.
    glycerol_set = set(fa_glycerol_carbons)
    glycerol_set.add(phospho_carbon)
    if len(glycerol_set) != 3:
        return False, "Fatty acid ester groups and the phosphate ester do not attach to three distinct glycerol carbons"
    
    # Now verify that these three carbons are adjacent (as expected in glycerol). 
    # In glycerol the two termini (the CH2 groups) are usually separated by the central CH.
    # We compute the pairwise graph distances between these three atoms and require they are ≤2.
    dmat = Chem.GetDistanceMatrix(mol)
    for i in glycerol_set:
        for j in glycerol_set:
            if i == j:
                continue
            if dmat[i][j] > 2:
                return False, "The three candidate glycerol carbons are not sufficiently connected (graph distance >2)"
    
    # -------------------------
    # Criterion 4: Check molecular weight (sanity check)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a phosphatidic acid"

    # If all criteria are met, we presume the molecule is a phosphatidic acid.
    return True, "Contains a phosphate group and exactly two fatty acid esters attached to a connected three-carbon glycerol backbone"

# Example usage:
if __name__ == "__main__":
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Is phosphatidic acid?", result)
    print("Reason:", reason)