"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nucleobase and a ribose or deoxyribose sugar linked by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ribose/deoxyribose sugar pattern (5-membered ring with 1 oxygen)
    sugar_pattern = Chem.MolFromSmarts("[CX4]1[OX2][CX4][CX4][CX4]1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No ribose or deoxyribose sugar found."
    
    # More specific check for carbon 2 and 3
    sugar_atoms = set()
    for match in sugar_matches:
        for a_idx in match:
            sugar_atoms.add(a_idx)
    
    found_sugar_valid = False
    for match in sugar_matches:
       c1 = match[0]
       c2 = match[2]
       c3 = match[3]
       
       c2_neigh = set([n.GetIdx() for n in mol.GetAtomWithIdx(c2).GetNeighbors()])
       c3_neigh = set([n.GetIdx() for n in mol.GetAtomWithIdx(c3).GetNeighbors()])

       if any (mol.GetAtomWithIdx(n).GetAtomicNum() == 8 for n in c2_neigh) and any(mol.GetAtomWithIdx(n).GetAtomicNum() == 8 for n in c3_neigh):
           found_sugar_valid = True
           break
       elif any(mol.GetAtomWithIdx(n).GetAtomicNum() == 8 for n in c3_neigh):
           found_sugar_valid = True
           break
    
    if not found_sugar_valid:
       return False, "The ribose or deoxyribose sugar is not correctly formed."


    # Nucleobase patterns (purine or pyrimidine) - more general patterns
    purine_pattern = Chem.MolFromSmarts("n1cncc2nc[nH]c12")
    pyrimidine_pattern = Chem.MolFromSmarts("n1ccnc(=O)[nH]1")

    purine_matches = mol.GetSubstructMatches(purine_pattern)
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_pattern)

    if not purine_matches and not pyrimidine_matches:
      return False, "No purine or pyrimidine nucleobase found."

    # Check for glycosidic bond (C-N bond linking sugar and base).
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[NX3]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_matches) != 1:
       return False, f"Found {len(glycosidic_matches)} C-N glycosidic bonds, must have exactly 1"
    
    #Check that the glycosidic bond carbon is part of the sugar
    
    connected = False
    for match in glycosidic_matches:
      glycosidic_c = match[0]
      if glycosidic_c in sugar_atoms:
        c1_neigh = set([n.GetIdx() for n in mol.GetAtomWithIdx(glycosidic_c).GetNeighbors()])
        if any (mol.GetAtomWithIdx(n).GetAtomicNum() == 8 for n in c1_neigh):
            connected = True
            break
    if not connected:
        return False, "Sugar and base not connected via glycosidic bond or incorrect connection"
    
    return True, "Contains a nucleobase linked to a ribose or deoxyribose sugar."