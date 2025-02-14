"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid contains both a carbohydrate and a lipid moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for carbohydrate moiety (pyranose/furanose ring with glycosidic linkages)
    pyranose_pattern = Chem.MolFromSmarts("[CX4]1[OX2][CX4][CX4][CX4][CX4]1")
    furanose_pattern = Chem.MolFromSmarts("[CX4]1[OX2][CX4][CX4][CX4]1")
    glycosidic_linkage = Chem.MolFromSmarts("[OX2][CX4]")

    has_carbohydrate = False
    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        if len(mol.GetSubstructMatches(glycosidic_linkage)) >= 1: # Require at least one glycosidic linkage
           has_carbohydrate = True

    if not has_carbohydrate:
      # No carbohydrate detected
        return False, "No carbohydrate moiety detected (no pyranose/furanose ring with glycosidic bond)"

    # 2. Check for lipid moiety (long hydrocarbon chain or lipid-associated functional groups connected to saccharide)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") # at least 8 C's in a row.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([OX1])([OX1])")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)[OX1]")

    has_lipid = False
    
    # Check for connection between carbohydrate and lipid
    saccharide_atoms = mol.GetSubstructMatches(pyranose_pattern) if mol.HasSubstructMatch(pyranose_pattern) else mol.GetSubstructMatches(furanose_pattern)
    if len(saccharide_atoms) > 0: # Check connection only if carbohydrate part exists
        saccharide_atoms = [atom_idx for match in saccharide_atoms for atom_idx in match]

        # Check if any lipid moiety is directly bonded to the saccharide.
        for atom_idx in saccharide_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
               if neighbor.GetSymbol() == 'C':
                 # Check if this carbon is part of long chain
                 if mol.HasSubstructMatch(long_chain_pattern, [neighbor.GetIdx()]):
                     has_lipid = True
                     break
               if neighbor.GetSymbol() == 'O':
                    for neighbor_of_o in neighbor.GetNeighbors():
                       if neighbor_of_o.GetSymbol() == 'C':
                         # Check if this carbon is part of ester, ether or amide etc
                         if mol.HasSubstructMatch(ester_pattern, [neighbor_of_o.GetIdx()]) or mol.HasSubstructMatch(amide_pattern, [neighbor_of_o.GetIdx()]) or mol.HasSubstructMatch(phosphate_pattern,[neighbor_of_o.GetIdx()]) or mol.HasSubstructMatch(sulfate_pattern, [neighbor_of_o.GetIdx()]):
                             has_lipid = True
                             break
            if has_lipid:
               break # Stop at first link

    if not has_lipid:
        return False, "No lipid moiety detected directly linked to the carbohydrate"

    # Check number of carbons - must be >10
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 12: # increased min carbons
        return False, "Too few carbons for a saccharolipid"

    # Check molecular weight - saccharolipids are typically large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450: # increased min molecular weight
       return False, "Molecular weight too low for saccharolipid"

    # 4. Return True and reason
    return True, "Contains both a carbohydrate and a lipid moiety"