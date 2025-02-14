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
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1[CH](O)[CH](O)[C@@H]([CH2X4,CHX4,CX4]O)O1"),  # ribose with correct stereochemistry
        Chem.MolFromSmarts("[C@H]1[CH](O)[CH](O)[C@@H]([CH2X4,CHX4,CX4])O1"),   # ribose with missing O
        Chem.MolFromSmarts("[C@H]1[CH](O)[CH](O)[C@@H](C)O1"),  # deoxyribose with correct stereochemistry
        Chem.MolFromSmarts("[C@H]1[CH](O)[CH](O)[C@@H]([CH3X4,CH2X4,CHX4,CX4])O1"), # general deoxy ribose
        Chem.MolFromSmarts("[C@@H]1[CH](O)[CH](O)[CH]([CH2X4,CHX4,CX4]O)O1"),  #  modified ribose with inverted stereocenter
         Chem.MolFromSmarts("[C@@H]1[CH](O)[CH](O)[CH]([CH2X4,CHX4,CX4])O1"),   # modified ribose with inverted stereocenter and missing O
         Chem.MolFromSmarts("[C@@H]1[CH](O)[CH](O)[CH](C)O1"),  #  modified deoxyribose with inverted stereocenter
        Chem.MolFromSmarts("[C@@H]1[CH](O)[CH](O)[CH]([CH3X4,CH2X4,CHX4,CX4])O1"),# general modified deoxy ribose
    ]

    sugar_matches = []
    for pattern in sugar_patterns:
      sugar_matches.extend(mol.GetSubstructMatches(pattern))


    if not sugar_matches:
        return False, "No sugar ring found"

    # 2. Check for a nucleobase (pyrimidine or purine)
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1[cn][cn][c](=[O,S])[nh]1"), # pyrimidines: modified or not.
        Chem.MolFromSmarts("n1[cn][cn][c][nh]1"), # modified pyrimidines
         Chem.MolFromSmarts("n1[cn][cn][c](=[O,S])[n]1"), # modified pyrimidines
        Chem.MolFromSmarts("n1cnc2[c]1[n][c][n][c]2[N,O]"), # purines: modified or not
        Chem.MolFromSmarts("n1cnc2[c]1[n][c][n][c]2"), # modified purines
         Chem.MolFromSmarts("n1cnc2[c]1[n][c][n][c]2[NH2,OH]"), #modified purines
    ]


    # 3. Check for phosphate
    phosphate_patterns = [
        Chem.MolFromSmarts("[OX2][P](=[OX1])([OX1,OX2])([OX1,OX2])"), #normal phosphate
       Chem.MolFromSmarts("[OX2]P([OX1])([OX1])O[CX4]"),  # phosphate connected to a carbon
        Chem.MolFromSmarts("O=P1([OX2])O[CX4][CX4]O1"),   #cyclic phosphate
        Chem.MolFromSmarts("[CX4][P](=[OX1])([OX1,OX2])([OX1,OX2])") # C-P phosphate
    ]

    # 4. Validate connectivity
    is_nucleoside_phosphate_found = False

    for sugar_match in sugar_matches:
        sugar_atoms = set(sugar_match)
        found_base = False
        found_phosphate = False


        for base_pattern in nucleobase_patterns:
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
                     found_base = True
                     break #stop search
             if found_base:
                break
        if not found_base:
            continue #check other sugar

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
                        found_phosphate = True
                        break

                if found_phosphate:
                    break
            if found_phosphate:
                break
        if found_phosphate:
            is_nucleoside_phosphate_found = True
            break #break if we found a valid nucleoside phosphate

    if not is_nucleoside_phosphate_found:
        return False, "Molecule is not a valid nucleoside phosphate."


    return True, "Molecule is a valid nucleoside phosphate."