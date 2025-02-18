"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: CHEBI:24676 nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate consists of a nucleobase, a sugar, and one or more phosphate groups,
    with specific connectivity.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for phosphate groups (at least one)
    # allow for mono, di, tri phosphates and also phosphonates/cyclic phosphates
    phosphate_pattern = Chem.MolFromSmarts("[OX1,OX2][P](=[OX1])([OX1,OX2])([OX1,OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
         return False, "No phosphate group found"

    # 2. Check for a sugar ring (ribose, deoxyribose, or modified versions)
    # More general pattern: 5-membered ring with at least one oxygen and 2 hydroxyls
    sugar_pattern = Chem.MolFromSmarts("C1[CH](O)[CH](O)[CH]([CH2X4,CHX4,CX4]O)O1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
      return False, "No sugar ring found"


    # 3. Check for a nucleobase (pyrimidine or purine)
    # Pyrimidines (cytosine, thymine, uracil)
    pyrimidine_pattern1 = Chem.MolFromSmarts("n1cncc(=O)[nH]1") #cytosine/uracil
    pyrimidine_pattern2 = Chem.MolFromSmarts("n1cc(C)c(=O)[nH]1") #thymine
    pyrimidine_pattern3 = Chem.MolFromSmarts("n1cncn1=O") # modified pyrimidine
    pyrimidine_pattern4 = Chem.MolFromSmarts("n1cncc[nH]1") #modified cytosine
    # Purines (adenine, guanine, hypoxanthine)
    purine_pattern1 = Chem.MolFromSmarts("n1cnc2c1ncnc2N") #adenine
    purine_pattern2 = Chem.MolFromSmarts("n1cnc2c1nc(=O)n[nH]2") #guanine
    purine_pattern3 = Chem.MolFromSmarts("n1cnc2c1ncnc2O") #hypoxanthine
    purine_pattern4 = Chem.MolFromSmarts("n1cncn2c1nc[nH]c2") # modified purine
    nucleobase_matches = mol.GetSubstructMatches(pyrimidine_pattern1)
    nucleobase_matches.extend(mol.GetSubstructMatches(pyrimidine_pattern2))
    nucleobase_matches.extend(mol.GetSubstructMatches(pyrimidine_pattern3))
    nucleobase_matches.extend(mol.GetSubstructMatches(pyrimidine_pattern4))
    nucleobase_matches.extend(mol.GetSubstructMatches(purine_pattern1))
    nucleobase_matches.extend(mol.GetSubstructMatches(purine_pattern2))
    nucleobase_matches.extend(mol.GetSubstructMatches(purine_pattern3))
    nucleobase_matches.extend(mol.GetSubstructMatches(purine_pattern4))
    if not nucleobase_matches:
        return False, "No nucleobase found"
    
    # 4. Check connectivity: the phosphate must be attached to the sugar, and the nucleobase to the sugar.
    # We'll use the atom indices found in substructure matches to verify the connectivity.
    
    is_connected = False
    for sugar_match in sugar_matches:
        for phosphate_match in phosphate_matches:
            for base_match in nucleobase_matches:
               #Check phosphate is connected to sugar
                phosphate_atoms = set(phosphate_match)
                sugar_atoms = set(sugar_match)
                found_connection = False
                for p_atom_idx in phosphate_atoms:
                  p_atom = mol.GetAtomWithIdx(p_atom_idx)
                  for neighbor in p_atom.GetNeighbors():
                    if neighbor.GetIdx() in sugar_atoms:
                        found_connection = True
                        break #phosphate connected to sugar
                  if found_connection:
                     break
                if not found_connection:
                  continue #check other phosphate group

                #Check if the nucleobase is connected to sugar

                base_atoms = set(base_match)
                found_base_sugar_connection = False

                for b_atom_idx in base_atoms:
                  b_atom = mol.GetAtomWithIdx(b_atom_idx)
                  for neighbor in b_atom.GetNeighbors():
                    if neighbor.GetIdx() in sugar_atoms:
                        found_base_sugar_connection = True
                        break # base connected to sugar
                  if found_base_sugar_connection:
                     break
                if not found_base_sugar_connection:
                    continue  # check next nucleobase match

                
                is_connected = True
                break #break if found one valid nucleoside phosphate
            if is_connected:
               break
        if is_connected:
            break


    if not is_connected:
        return False, "No valid nucleoside phosphate connectivity"

    return True, "Contains a nucleobase, a sugar, and at least one phosphate group with correct connectivity"