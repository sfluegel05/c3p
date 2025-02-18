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

    # 1. Check for sugar ring (ribose, deoxyribose, or modified versions)
    # More general pattern: 5-membered ring with at least one oxygen and 2 hydroxyls.
    # Note that this can be modified
    sugar_pattern = Chem.MolFromSmarts("C1[CH](O)[CH](O)[CH]([CH2X4,CHX4,CX4]O)O1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar ring found"


    # 2. Check for a nucleobase (pyrimidine or purine)
    # Pyrimidines (cytosine, thymine, uracil)
    pyrimidine_patterns = [
        Chem.MolFromSmarts("n1cncc(=O)[nH]1"),  # cytosine/uracil
        Chem.MolFromSmarts("n1cc(C)c(=O)[nH]1"),  # thymine
        Chem.MolFromSmarts("n1cncn1=O"),  # modified pyrimidine
        Chem.MolFromSmarts("n1cncc[nH]1"), #modified cytosine
    ]
    # Purines (adenine, guanine, hypoxanthine)
    purine_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2N"),  # adenine
        Chem.MolFromSmarts("n1cnc2c1nc(=O)n[nH]2"),  # guanine
        Chem.MolFromSmarts("n1cnc2c1ncnc2O"), #hypoxanthine
         Chem.MolFromSmarts("n1cncn2c1nc[nH]c2"), # modified purine
    ]

    #3. Check for phosphate
    phosphate_patterns = [
        Chem.MolFromSmarts("[OX2][P](=[OX1])([OX1,OX2])([OX1,OX2])"), #normal phosphate
        Chem.MolFromSmarts("[CX4][P](=[OX1])([OX1,OX2])([OX1,OX2])") # C-P phosphate
    ]

    #4. Validate connectivity
    is_nucleoside_phosphate_found = False
    for sugar_match in sugar_matches:
        sugar_atoms = set(sugar_match)

        for base_pattern in pyrimidine_patterns + purine_patterns:
             base_matches = mol.GetSubstructMatches(base_pattern)
             if not base_matches:
                 continue
             for base_match in base_matches:
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
                     continue #check other base pattern
                 for phosphate_pattern in phosphate_patterns:
                     phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
                     if not phosphate_matches:
                        continue
                     for phosphate_match in phosphate_matches:
                         phosphate_atoms = set(phosphate_match)
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
                         is_nucleoside_phosphate_found=True
                         break # break if found one valid nucleoside phosphate
                     if is_nucleoside_phosphate_found:
                       break
                 if is_nucleoside_phosphate_found:
                    break
             if is_nucleoside_phosphate_found:
                break
        if is_nucleoside_phosphate_found:
           break


    if not is_nucleoside_phosphate_found:
        return False, "Molecule is not a valid nucleoside phosphate."

    return True, "Molecule is a valid nucleoside phosphate."