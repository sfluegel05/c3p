"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: Phosphatidic acid (PA)
Definition:
  A derivative of glycerol in which one hydroxy group (commonly but not necessarily primary) is esterified with phosphoric acid 
  and the other two are esterified with fatty acids.
Improvements in this version:
  - Identify a phosphate group by checking for a P atom with ≥3 oxygen neighbors and exactly one oxygen that attaches to a carbon (the phosphate ester).
  - Exclude phosphate groups that carry any explicit negative formal charge on the P or one of its immediate oxygen neighbors.
  - Identify fatty acid ester groups with the substructure "OC(=O)" while filtering out any bridging oxygens attached directly to phosphorus.
  - For each fatty acid ester find the glycerol carbon (avoiding the carbonyl carbon with a C=O) and then try to identify the phosphate-derived glycerol carbon.
  - Verify that the three candidate glycerol carbons (2 from FA esters and 1 from the phosphate ester) are distinct and contiguous (graph distance ≤2).
  - Do a basic molecular weight sanity check.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid (PA) based on its SMILES string.
    PA is a glycerol derivative where one hydroxy is esterified with phosphoric acid and the other two with fatty acids.
    
    Args:
      smiles (str): Input SMILES for the molecule.
    
    Returns:
      bool: True if the molecule qualifies as a phosphatidic acid, otherwise False.
      str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------------------------
    # Criterion 1: Identify a suitable phosphate group.
    # Look for a phosphorus (P) atom that has at least three oxygen neighbors.
    # Among those oxygens, exactly one should connect to a carbon (the phosphate–glycerol bridging).
    # Also, none of the atoms in the phosphate head (P and its O neighbors) should carry a negative formal charge.
    valid_P_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "P":
            continue
        # Check formal charge on the phosphorus.
        if atom.GetFormalCharge() < 0:
            continue
        
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"]
        if len(oxy_neighbors) < 3:
            continue
        
        # Eliminate any phosphate group with any negatively charged oxygen neighbor.
        neg_charge = False
        for o in oxy_neighbors:
            if o.GetFormalCharge() < 0:
                neg_charge = True
                break
        if neg_charge:
            continue
        
        # Count among these oxygens how many are attached to at least one carbon (other than the P).
        o_attached_to_C = 0
        for o in oxy_neighbors:
            for nbr2 in o.GetNeighbors():
                if nbr2.GetIdx() == atom.GetIdx():
                    continue
                if nbr2.GetSymbol() == "C":
                    o_attached_to_C += 1
                    break
        if o_attached_to_C == 1:
            valid_P_atom = atom
            break
    if valid_P_atom is None:
        return False, "No phosphate group found with ≥3 oxygen neighbors (all neutral) and exactly one oxygen attached to a carbon"
    
    # -------------------------
    # Criterion 2: Identify fatty acid ester (FA) groups.
    # Look for the substructure "OC(=O)". Then make sure that the bridging oxygen is not directly attached to phosphorus.
    fa_pattern = Chem.MolFromSmarts("OC(=O)")
    fa_matches = mol.GetSubstructMatches(fa_pattern)
    fa_candidates = []    # To collect candidate glycerol carbon indices (from the FA ester linkage)
    bridging_O_indices = set()  # To mark oxygens used in FA linkages.
    
    for match in fa_matches:
        # Each match is a tuple: (bridging O, carbonyl C, carbonyl O)
        bridging_O_idx = match[0]
        bridging_O = mol.GetAtomWithIdx(bridging_O_idx)
        # Skip if bridging oxygen is bonded to any phosphorus.
        if any(nbr.GetSymbol() == "P" for nbr in bridging_O.GetNeighbors()):
            continue
        
        # From the bridging O, choose a neighbor that is a carbon and is not the carbonyl carbon.
        candidate = None
        for nbr in bridging_O.GetNeighbors():
            if nbr.GetSymbol() != "C":
                continue
            # Distinguish the glycerol carbon from the carbonyl carbon by checking for a double bonded oxygen.
            has_double_bond_oxygen = False
            for nbr2 in nbr.GetNeighbors():
                if nbr2.GetIdx() == bridging_O.GetIdx():
                    continue
                bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                if nbr2.GetSymbol() == "O" and bond is not None and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_double_bond_oxygen = True
                    break
            if not has_double_bond_oxygen:
                candidate = nbr.GetIdx()
                break
        if candidate is not None:
            fa_candidates.append(candidate)
            bridging_O_indices.add(bridging_O_idx)
    if len(fa_candidates) != 2:
        return False, f"Found {len(fa_candidates)} fatty acid ester group(s) after filtering; expected exactly 2."
    
    # -------------------------
    # Criterion 3: Identify the glycerol carbon from the phosphate ester.
    # Look at the P atom’s oxygen neighbors (excluding those already used for FA ester linkages)
    # and find candidate oxygens that lead to a carbon.
    phospho_candidate_list = []
    for o in valid_P_atom.GetNeighbors():
        if o.GetSymbol() != "O":
            continue
        if o.GetIdx() in bridging_O_indices:
            continue
        # For each such oxygen, if it leads to a carbon (other than the P) then add that carbon as candidate.
        for nbr in o.GetNeighbors():
            if nbr.GetIdx() == valid_P_atom.GetIdx():
                continue
            if nbr.GetSymbol() == "C":
                phospho_candidate_list.append(nbr.GetIdx())
                break   # Taking the first carbon is sufficient for this oxygen.
                
    # We expect the phosphate ester to contribute exactly one glycerol carbon.
    # But if more than one candidate is found, test each one by forming a set with the two FA candidates.
    candidate_found = None
    dmat = Chem.GetDistanceMatrix(mol)
    for pc in phospho_candidate_list:
        glycerol_set = set(fa_candidates)
        glycerol_set.add(pc)
        if len(glycerol_set) != 3:
            continue
        # Check connectivity: graph distance between any two glycerol candidates must be ≤ 2.
        distances_ok = True
        for a in glycerol_set:
            for b in glycerol_set:
                if a == b:
                    continue
                if dmat[a][b] > 2:
                    distances_ok = False
                    break
            if not distances_ok:
                break
        if distances_ok:
            candidate_found = pc
            break
    if candidate_found is None:
        return False, "Could not find a single phosphate-derived glycerol carbon that connects with the two FA ester glycerol carbons within a graph distance of 2."
    
    # Now the three glycerol carbons are the two from FA ester groups and candidate_found.
    glycerol_set = set(fa_candidates)
    glycerol_set.add(candidate_found)
    if len(glycerol_set) != 3:
        return False, "Fatty acid ester groups and the phosphate ester do not attach to three distinct glycerol carbons."
    
    # -------------------------
    # Criterion 4: (Optional) Molecular weight sanity check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a phosphatidic acid"
    
    return True, "Contains a phosphate group (with neutral head) and exactly two fatty acid esters attached to a connected three‐carbon glycerol backbone"

# Example usage:
if __name__ == "__main__":
    # Test one of the provided true positive SMILES (phosphatidic acid with two FA esters and a phosphate ester).
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Is phosphatidic acid?", result)
    print("Reason:", reason)