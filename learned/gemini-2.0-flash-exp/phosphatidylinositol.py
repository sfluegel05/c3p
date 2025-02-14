"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    A phosphatidylinositol has a myo-inositol ring, a glycerol backbone with two fatty acids, and a phosphate linking the inositol and glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a phosphatidylinositol, False otherwise
         str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for Inositol Ring
    inositol_pattern = Chem.MolFromSmarts("[C]1([C]([C]([C]([C]([C]1O)O)O)O)O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # 2. Check for Glycerol Backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
       return False, "No glycerol backbone found"

    # 3. Check for Phosphate Group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)[OX1]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate group found"

    # 4. Verify that inositol and glycerol are connected via phosphate
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)

    connected = False
    for i_match in inositol_matches:
        i_atom = mol.GetAtomWithIdx(i_match[0]) # Inositol carbon
        for p_match in phosphate_matches:
            p_atom = mol.GetAtomWithIdx(p_match[0]) # Phosphate atom
            for g_match in glycerol_matches:
                for gly_atom_idx in g_match: #check each of the glycerol carbons
                  gly_atom = mol.GetAtomWithIdx(gly_atom_idx)
                  for nbr in gly_atom.GetNeighbors():
                      if nbr.GetAtomicNum() == 8 : # check that the glycerol oxygen is connected to phosphate
                           if mol.GetBondBetweenAtoms(nbr.GetIdx(), p_atom.GetIdx()):
                                # check if inositol ring is connected to phosphate oxygen
                                for nbr_p in p_atom.GetNeighbors():
                                    if nbr_p.GetAtomicNum() == 8:
                                        for nbr_i in i_atom.GetNeighbors():
                                            if nbr_i.GetAtomicNum() == 8:
                                                if mol.GetBondBetweenAtoms(nbr_p.GetIdx(), nbr_i.GetIdx()):
                                                    connected=True
                                                    break
                                        if connected:
                                          break
                  if connected:
                    break

                if connected:
                   break
            if connected:
                break
        if connected:
            break

    if not connected:
        return False, "Inositol and glycerol are not connected via phosphate"

    # 5. Check for two fatty acid chains attached via esters.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
      return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    if len(fatty_acid_matches) < 2 :
       return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"


    # Check for number of rotatable bonds to exclude short chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8 :
        return False, "Chains too short to be fatty acids"


    return True, "Contains inositol, glycerol, phosphate, and two fatty acid chains"