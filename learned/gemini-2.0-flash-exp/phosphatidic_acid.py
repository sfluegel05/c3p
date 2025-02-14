"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid has a glycerol backbone, a phosphate group, and two fatty acid chains.
    The two fatty acid chains are attached via ester bonds to the glycerol.
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
         return False, "No glycerol backbone found"
    
    # Find phosphate group
    phosphate_pattern = Chem.MolFromSmarts("[P](=[OX1])([OX2])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
         return False, "No phosphate group found"
    
    # Verify attachment between glycerol and phosphate
    glycerol_oxygen_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-O")
    attached = False
    for glycerol_match in glycerol_matches:
      for atom_idx in glycerol_match:
        glycerol_atom = mol.GetAtomWithIdx(atom_idx)
        if glycerol_atom.GetAtomicNum() == 6: #verify carbon is from glycerol
          for neighbor in glycerol_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
               for phosphate_match in phosphate_matches:
                 for p_idx in phosphate_match:
                   phosphate_atom = mol.GetAtomWithIdx(p_idx)
                   if phosphate_atom.GetAtomicNum() == 15:
                     for p_neighbor in phosphate_atom.GetNeighbors():
                       if p_neighbor.GetIdx() == neighbor.GetIdx(): #verify that phosphate is connected by an oxygen
                           attached = True
                           break
                 if attached:
                    break
          if attached:
             break
      if attached:
          break

    if not attached:
      return False, "Phosphate not directly attached to glycerol"

    # Look for ester groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Verify exactly two ester groups are connected to glycerol oxygens
    ester_count = 0
    for glycerol_match in glycerol_matches:
        for atom_idx in glycerol_match:
            glycerol_atom = mol.GetAtomWithIdx(atom_idx)
            if glycerol_atom.GetAtomicNum() == 6: #verify carbon is from glycerol
               for neighbor in glycerol_atom.GetNeighbors(): #check the neighbors (oxygens) of glycerol carbons
                 if neighbor.GetAtomicNum() == 8:
                   for ester_match in ester_matches:
                        for ester_atom in ester_match:
                             ester_atom_obj = mol.GetAtomWithIdx(ester_atom)
                             if ester_atom_obj.GetAtomicNum() == 8: #check for ester oxygen
                                 if ester_atom_obj.GetIdx() == neighbor.GetIdx():
                                       ester_count += 1
    if ester_count != 2:
      return False, f"Found {ester_count} ester groups attached to glycerol, need exactly 2"

    # Check for fatty acid chains connected to the ester groups: search for >= 2 carbons
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3,CX4]~[CX3,CX4]")
    fatty_acid_count = 0
    for ester_match in ester_matches:
        for ester_atom in ester_match:
          ester_atom_obj = mol.GetAtomWithIdx(ester_atom)
          if ester_atom_obj.GetAtomicNum() == 6:
             for neighbor in ester_atom_obj.GetNeighbors():
                  if neighbor.GetAtomicNum() == 8:
                     for next_neighbor in ester_atom_obj.GetNeighbors():
                        if next_neighbor.GetAtomicNum() == 6:
                           
                            chain_match = mol.GetSubstructMatches(fatty_acid_pattern, [next_neighbor.GetIdx()])

                            if chain_match:
                                fatty_acid_count += 1
                                break
    
    if fatty_acid_count != 2:
       return False, f"Missing two fatty acid chains attached to the ester bonds, found {fatty_acid_count}"
   

    return True, "Contains glycerol backbone with a phosphate group and two fatty acid chains attached via ester bonds"