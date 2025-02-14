"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    These molecules have a glycerol backbone esterified at the sn-1 position with a fatty acid,
    and at the sn-3 position with a phosphoethanolamine group and (R) configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Glycerol Backbone with (R)-configuration at central carbon
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)(CO)CO") # [C@H](O)(C)-C-O
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
      return False, "No glycerol backbone with correct stereochemistry found"
    if len(glycerol_matches) != 1:
        return False, "More than one glycerol group found"

    glycerol_match = glycerol_matches[0]
    glycerol_ch2_1_idx = glycerol_match[2]  #sn-1
    glycerol_ch2_2_idx = glycerol_match[3] # sn-3


    # 2. Phosphoethanolamine Group
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    phosphoethanolamine_matches = mol.GetSubstructMatches(phosphoethanolamine_pattern)
    if not phosphoethanolamine_matches or len(phosphoethanolamine_matches) != 1:
       return False, "No phosphoethanolamine group found or found more than one"
    
    phosphoethanolamine_match = phosphoethanolamine_matches[0]
    phosphorus_idx = phosphoethanolamine_match[1]
    
    phospho_oxygen_glycerol_idx = -1
    for neighbor in mol.GetAtomWithIdx(phosphorus_idx).GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            for neighbor2 in mol.GetAtomWithIdx(neighbor.GetIdx()).GetNeighbors():
                if neighbor2.GetIdx() == glycerol_ch2_2_idx:
                    phospho_oxygen_glycerol_idx = neighbor.GetIdx()
                    break
            if phospho_oxygen_glycerol_idx != -1:
                break
    if phospho_oxygen_glycerol_idx == -1:
       return False, "Phosphate group not correctly attached to the glycerol backbone at sn-3"

    # 3. Ester at sn-1
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    ester_at_sn1 = False
    if ester_matches:
        for match in ester_matches:
            carbonyl_carbon_idx = match[0]
            ester_oxygen_idx = match[1]
            
            for neighbor in mol.GetAtomWithIdx(ester_oxygen_idx).GetNeighbors():
               if neighbor.GetIdx() == glycerol_ch2_1_idx:
                   ester_at_sn1 = True
                   break
            if ester_at_sn1:
              break
    if not ester_at_sn1:
        return False, "No ester group at sn-1 glycerol position"
    
    # 4. Fatty Acid Chain (checking a minimum length of 3 carbons)
    fatty_acid_connected = False
    if ester_matches:
        for match in ester_matches:
          carbonyl_carbon_idx = match[0]
          for neighbor in mol.GetAtomWithIdx(carbonyl_carbon_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                chain_len = 0
                current_atom_idx = neighbor.GetIdx()
                visited_atoms = {carbonyl_carbon_idx, ester_oxygen_idx}
                while True:
                    
                  chain_len += 1
                  neighbors = [n for n in mol.GetAtomWithIdx(current_atom_idx).GetNeighbors() if n.GetIdx() not in visited_atoms and n.GetAtomicNum() == 6]
                  if len(neighbors) != 1:
                      break
                  
                  visited_atoms.add(current_atom_idx)
                  current_atom_idx = neighbors[0].GetIdx()

                if chain_len >= 3:
                    fatty_acid_connected = True
                    break
        if not fatty_acid_connected:
             return False, "Missing or too short fatty acid chain"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
         return False, "Molecular weight too low for a phosphoethanolamine"

    return True, "Matches 1-acyl-sn-glycero-3-phosphoethanolamine criteria"