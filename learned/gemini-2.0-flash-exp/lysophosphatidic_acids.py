"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a glycerol backbone with one fatty acid chain attached via an ester bond and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens, where one will be connected to phosphate and one to the acyl chain)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for phosphate group (-P(=O)(O)-OH or -P(=O)(O)2
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([OX1])([OX1])[OX2]") #Modified to include a potentially missing oxygen
    if not mol.HasSubstructMatch(phosphate_pattern):
       return False, "No phosphate group found"
    
    # Check that the glycerol backbone is linked to the phosphate group
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    phosphate_match = mol.GetSubstructMatch(phosphate_pattern)

    if not glycerol_match or not phosphate_match:
        return False, "Could not find both glycerol and phosphate"


    #Look for a bond between the glycerol and the phosphate
    found_link = False
    for g_atom_idx in glycerol_match:
        g_atom = mol.GetAtomWithIdx(g_atom_idx)
        for p_atom_idx in phosphate_match:
            p_atom = mol.GetAtomWithIdx(p_atom_idx)
            for neighbor in g_atom.GetNeighbors():
                 if neighbor.GetIdx() == p_atom_idx:
                      found_link = True
                      break
        if found_link:
            break
    if not found_link:
       return False, "Phosphate group not connected to glycerol backbone"

    # Look for 1 ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Look for a fatty acid chain (long carbon chain attached to the ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
         return False, f"No fatty acid chain found"
    
    #Check that fatty acid chain is connected to the ester group
    found_ester_link = False
    for e_match in ester_matches:
         ester_atom_idx = e_match[0]
         ester_atom = mol.GetAtomWithIdx(ester_atom_idx)
         for fa_match in fatty_acid_matches:
              for fa_atom_idx in fa_match:
                  fa_atom = mol.GetAtomWithIdx(fa_atom_idx)
                  for neighbor in ester_atom.GetNeighbors():
                      if neighbor.GetIdx() == fa_atom_idx:
                          found_ester_link = True
                          break
              if found_ester_link:
                  break

         if found_ester_link:
             break
    
    if not found_ester_link:
         return False, "Fatty acid chain not connected to ester group"

    # Check the number of phosphorus atoms - must be 1.
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
         return False, f"Must have exactly 1 phosphorus atom, found {p_count}"
    
    # Count rotatable bonds to verify long chains.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
       return False, "Chain too short to be a fatty acid"
    
    # Check molecular weight - LPA typically > 300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a lysophosphatidic acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 8:
         return False, "Too few carbons for a lysophosphatidic acid"
    if o_count < 6:
         return False, "Must have at least 6 oxygens"

    return True, "Contains a glycerol backbone with a single fatty acid chain and a phosphate group"