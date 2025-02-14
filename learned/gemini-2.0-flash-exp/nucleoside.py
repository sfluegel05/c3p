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
    
    # Nucleobase patterns (purine or pyrimidine)
    purine_pattern = Chem.MolFromSmarts("[n]1cncn2c1[n]c[n]2")
    pyrimidine_pattern = Chem.MolFromSmarts("[n]1cc[n]c(=[OX1])c1")

    purine_matches = mol.GetSubstructMatches(purine_pattern)
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_pattern)

    if not purine_matches and not pyrimidine_matches:
        return False, "No purine or pyrimidine nucleobase found."

    # Check for glycosidic bond (C-N bond linking sugar and base)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[NX3]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_matches:
        return False, "No C-N glycosidic bond found between sugar and base"
    
    # Check that the sugar and base are actually connected to each other via glycosidic bond
    sugar_atoms = set()
    for match in sugar_matches:
        for a_idx in match:
            sugar_atoms.add(a_idx)

    base_atoms = set()
    if purine_matches:
        for match in purine_matches:
            for a_idx in match:
                base_atoms.add(a_idx)
    if pyrimidine_matches:
        for match in pyrimidine_matches:
            for a_idx in match:
                base_atoms.add(a_idx)

    glycosidic_atoms = set()
    for match in glycosidic_matches:
        for a_idx in match:
            glycosidic_atoms.add(a_idx)

    connected = False
    for a_idx in glycosidic_atoms:
        neighbor_atoms = [neighbor.GetIdx() for neighbor in mol.GetAtomWithIdx(a_idx).GetNeighbors()]
        for neighbor in neighbor_atoms:
            if neighbor in sugar_atoms:
              for neighbor2 in neighbor_atoms:
                if neighbor2 in base_atoms:
                  connected=True
                  break
            if connected:
              break
        if connected:
            break

    if not connected:
      return False, "Sugar and base not connected via glycosidic bond"
    
    # Additional check on number of carbons and oxygen in the sugar
    sugar_c_count = 0
    sugar_o_count = 0
    for a_idx in sugar_atoms:
        atom = mol.GetAtomWithIdx(a_idx)
        if atom.GetAtomicNum() == 6:
            sugar_c_count += 1
        elif atom.GetAtomicNum() == 8:
            sugar_o_count += 1
    if not (4 <= sugar_c_count <= 5 and 3 <= sugar_o_count <= 4):
        return False, "Incorrect carbon or oxygen count on sugar"
    
    return True, "Contains a nucleobase linked to a ribose or deoxyribose sugar."