"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide is a phosphatidylinositol that is phosphorylated at one or more of the hydroxy groups of inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with a phosphate group attached via an oxygen
    # Using a slightly more relaxed SMARTS pattern to match both P=O and P-O groups
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2][CX3](=[OX1,~0]))[CH2X4][OX2][P](=[OX1,~0])([OX2])([OX2])")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
         return False, "No glycerol-phosphate backbone found"

    # Check for a myo-inositol ring (six membered ring with 5 OH groups)
    inositol_pattern = Chem.MolFromSmarts("C1[C]([O])[C]([O])[C]([O])[C]([O])[C]([O])1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol ring found"

    # Check for a phosphate group on the inositol ring, directly or indirectly linked via an oxygen
    # Modified to allow for linkers between inositol and phosphate
    inositol_phosphate_pattern = Chem.MolFromSmarts("[C]~[O]~[P](=[OX1,~0])([OX2])([OX2])")
    inositol_phosphate_matches = mol.GetSubstructMatches(inositol_phosphate_pattern)

    # Exclude the glycerol phosphate by finding the match linked to glycerol
    glycerol_phosphate_match = mol.GetSubstructMatch(glycerol_phosphate_pattern)
    if glycerol_phosphate_match:
      glycerol_phosphate_atoms = [glycerol_phosphate_match[8]] # the P atom of the glycerol phosphate
    else:
      glycerol_phosphate_atoms = []
    
    # Check that a phosphate group in the match is linked to the inositol ring and is not the same as the one on glycerol
    additional_phosphate_count = 0
    for match in inositol_phosphate_matches:
      phosphorus_atom_idx = match[2] # the P atom in the inositol linked phosphate group.
      phosphorus_atom = mol.GetAtomWithIdx(phosphorus_atom_idx)
      
      # Find the inositol carbon connected to the phosphate
      found_inositol_connection = False
      for neighbor in phosphorus_atom.GetNeighbors():
          for inositol_atom in mol.GetSubstructMatch(inositol_pattern):
            if neighbor.GetIdx() == inositol_atom:
                found_inositol_connection = True
                break
          if found_inositol_connection:
            break

      # Only count this phosphate if it is not the one on glycerol and is linked to inositol
      if found_inositol_connection and (not (phosphorus_atom_idx in glycerol_phosphate_atoms)):
        additional_phosphate_count += 1

    if additional_phosphate_count < 1:
        return False, "Less than one additional phosphate group found on the inositol ring"

    return True, "Contains glycerol backbone, a phosphate group linked to glycerol, a myo-inositol ring, and at least one additional phosphate group on the inositol ring"