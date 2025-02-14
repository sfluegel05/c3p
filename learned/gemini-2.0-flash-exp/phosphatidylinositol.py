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

    # 2. Check for Glycerol Backbone with phosphate group
    glycerol_phosphate_pattern = Chem.MolFromSmarts("CC(COP(=O)(O)[O])C")
    glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if len(glycerol_phosphate_matches) == 0:
        return False, "No glycerol backbone with phosphate group found"

    # 3. Verify that inositol is connected to glycerol via the phosphate
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)

    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    
    connected = False
    for p_match in phosphate_matches:
        for i_match in inositol_matches:
            for glycerol_match in glycerol_phosphate_matches:
                p_atom = mol.GetAtomWithIdx(p_match[0])
                i_atom = mol.GetAtomWithIdx(i_match[0])
                gly_atom_list = [mol.GetAtomWithIdx(glycerol_match[0]), mol.GetAtomWithIdx(glycerol_match[1]), mol.GetAtomWithIdx(glycerol_match[2]), mol.GetAtomWithIdx(glycerol_match[3])] #list of atoms in glycerol
                
                for gly_atom in gly_atom_list: #check if each atom of glycerol has bond to phosphate
                  if mol.GetBondBetweenAtoms(p_atom.GetIdx(), gly_atom.GetIdx()):
                      for nbr in p_atom.GetNeighbors():  #check if phosphate has bond to inositol
                          if nbr.GetIdx() == i_atom.GetIdx():
                             connected = True
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
       return False, "Inositol is not connected to glycerol via phosphate"

    # 4. Check for two fatty acid chains attached via esters.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
      return False, f"Found {len(ester_matches)} ester groups, need at least 2"


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