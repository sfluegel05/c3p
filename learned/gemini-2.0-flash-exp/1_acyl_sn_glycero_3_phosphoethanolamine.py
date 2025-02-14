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

    # 1. Glycerol Backbone
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)(CO)")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "No glycerol backbone with correct stereochemistry found"
    glycerol_carbon_idx = glycerol_match[0]
    glycerol_ch2_1_idx = -1
    glycerol_ch2_2_idx = -1

    for neighbor in mol.GetAtomWithIdx(glycerol_carbon_idx).GetNeighbors():
         if neighbor.GetAtomicNum() == 6 and len(mol.GetAtomWithIdx(neighbor.GetIdx()).GetNeighbors()) == 3:
            for neighbor2 in mol.GetAtomWithIdx(neighbor.GetIdx()).GetNeighbors():
                if neighbor2.GetAtomicNum() == 8:
                     glycerol_ch2_1_idx = neighbor.GetIdx()
                     break
         if neighbor.GetAtomicNum() == 6 and len(mol.GetAtomWithIdx(neighbor.GetIdx()).GetNeighbors()) == 4:
            glycerol_ch2_2_idx = neighbor.GetIdx()

    if glycerol_ch2_1_idx == -1 or glycerol_ch2_2_idx == -1:
        return False, "Glycerol chain not correctly formed"

    # 2. Ester at sn-1
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

    # 3. Phosphoethanolamine Group
    phosphoethanolamine_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    phosphoethanolamine_match = mol.GetSubstructMatch(phosphoethanolamine_pattern)
    if not phosphoethanolamine_match:
       return False, "No phosphoethanolamine group found"

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
        return False, "Phosphate group not correctly attached to the glycerol backbone"


    # 4. Fatty Acid Chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
         return False, "Missing fatty acid chain"

    # 5. Verify that there is a carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No carbonyl group found"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
         return False, "Molecular weight too low for a phosphoethanolamine"

    # 6. Check for exactly one glycerol, ester, and phosphoethanolamine group
    if not glycerol_match:
        return False, "Incorrect number of glycerol groups"
    if not ester_at_sn1:
        return False, "Incorrect number of ester groups"
    if not phosphoethanolamine_match:
        return False, "Incorrect number of phosphoethanolamine groups"


    return True, "Matches 1-acyl-sn-glycero-3-phosphoethanolamine criteria"