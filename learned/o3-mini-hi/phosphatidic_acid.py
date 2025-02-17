"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: Phosphatidic acid
Definition: A derivative of glycerol in which one hydroxy group is esterified with phosphoric acid
and the other two hydroxy groups are esterified with fatty acids.
Improvements:
  - Only count fatty acid “OC(=O)” matches whose oxygen is not directly attached to phosphorus.
  - For a P atom to qualify as the PA phosphate, it must be bonded to at least three oxygens,
    of which exactly one is attached to a carbon (to avoid headgroup substitutions seen in PG/PE).
  - Finally, we also ensure that the two filtered fatty acid esters attach to two different glycerol carbons.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    Phosphatidic acid (PA) is defined here as a molecule that contains:
      1. A phosphate group where the phosphorus is bonded to at least three O atoms,
         and exactly one of these O atoms is linked to a carbon (i.e. the glycerol linkage).
      2. Exactly two fatty acid ester bonds (OC(=O)) that are bridges from the glycerol to fatty acids.
         (Extra OC(=O) groups in a fatty acid chain or attached to the phosphate are not counted.)
      3. The two fatty acid ester bonds must be attached to two distinct carbon atoms (implying
         a glycerol backbone with two esterified hydroxyl groups).
      4. The molecule should have a reasonable molecular weight (here, >300 Da).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is phosphatidic acid; False otherwise.
        str: Explanation string for the classification decision.
    """
    # Parse the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Criterion 1: Find a suitable phosphate group.
    # We iterate over all phosphorus atoms and look for one that has at least three oxygen neighbors,
    # and exactly one such neighbor is bonded to a carbon (the glycerol attachment).
    phosphate_ok = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "P":
            oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "O"]
            if len(oxy_neighbors) < 3:
                continue
            # Count how many O neighbors are attached to a carbon.
            o_attached_to_C = 0
            for o in oxy_neighbors:
                # For each oxygen, if it is bonded to any carbon (other than P) count it.
                for nbr in o.GetNeighbors():
                    if nbr.GetIdx() == atom.GetIdx():
                        continue
                    if nbr.GetSymbol() == "C":
                        o_attached_to_C += 1
                        break
            if o_attached_to_C == 1:
                phosphate_ok = True
                break
    if not phosphate_ok:
        return False, "Phosphate group not found with required substitution (should have >=3 O neighbors and exactly one O linked to C)"
    
    # Criterion 2: Count fatty acid ester groups.
    # Find matches for the ester substructure "OC(=O)". This pattern identifies an oxygen 
    # bonded to a carbonyl carbon. However, we want the bridging oxygen that comes from the glycerol.
    fatty_acid_pattern = Chem.MolFromSmarts("OC(=O)")
    fa_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    
    # Filter out matches where the bridging oxygen (first atom in the match) is directly attached to phosphorus.
    filtered_matches = []
    # Also record the glycerol carbon attached to the bridging O for a later connectivity check.
    glycerol_carbons = set()
    
    for match in fa_matches:
        bridging_O_idx = match[0]  # the oxygen in the pattern "O C(=O)"
        bridging_O = mol.GetAtomWithIdx(bridging_O_idx)
        # Exclude if this oxygen is attached directly to a phosphorus (it would then be part of the phosphate head)
        if any(neighbor.GetSymbol() == "P" for neighbor in bridging_O.GetNeighbors()):
            continue
        filtered_matches.append(match)
        # Identify the glycerol carbon the bridging oxygen is attached to.
        for nbr in bridging_O.GetNeighbors():
            if nbr.GetSymbol() == "C":
                glycerol_carbons.add(nbr.GetIdx())
    
    if len(filtered_matches) != 2:
        return False, f"Found {len(filtered_matches)} fatty acid ester group(s); exactly 2 are expected (after filtering out phosphate-linked esters)"
    
    # Criterion 3: Check that the two fatty acid bridges come from two different glycerol-carbons.
    if len(glycerol_carbons) != 2:
        return False, "Fatty acid ester groups do not arise from two distinct glycerol hydroxyl groups"
    
    # Criterion 4: Check overall molecular weight as a sanity measure.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low for phosphatidic acid ({mol_wt:.1f} Da)"
    
    return True, "Contains appropriate phosphate group and exactly two fatty acid ester chains attached to distinct glycerol hydroxyls"

# Example usage (this may be removed or modified as needed):
if __name__ == "__main__":
    test_smiles = "P(OC[C@H](OC(=O)CCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Is phosphatidic acid?", result)
    print("Reason:", reason)