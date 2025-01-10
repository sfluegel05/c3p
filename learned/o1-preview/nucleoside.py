"""
Classifies: CHEBI:33838 nucleoside
"""
"""
Classifies: nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nucleobase attached to a ribose or deoxyribose sugar via an N-glycosidic bond.

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

    # Define general SMARTS patterns for purine and pyrimidine nucleobases
    purine_smarts = 'c1nc2ncnc-2n1'  # Purine ring system
    pyrimidine_smarts = 'c1ccncn1'    # Pyrimidine ring system

    purine_pattern = Chem.MolFromSmarts(purine_smarts)
    pyrimidine_pattern = Chem.MolFromSmarts(pyrimidine_smarts)

    # Define a flexible SMARTS pattern for furanose sugars (ribose and deoxyribose)
    sugar_smarts = '[C@@H]1([O])[C@H]([O])[C@@H]([O])[C@H](CO)[O]1'  # Five-membered sugar ring
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)

    # Check for nucleobase presence
    has_purine = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_pattern)
    if not (has_purine or has_pyrimidine):
        return False, "No purine or pyrimidine nucleobase found"

    # Check for sugar moiety presence
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar moiety found"

    # Identify nucleobase nitrogen that forms N-glycosidic bond
    nucleobase_nitrogen_idx = None
    if has_purine:
        # In purines, N9 forms the glycosidic bond
        purine_matches = mol.GetSubstructMatches(purine_pattern)
        for match in purine_matches:
            n9_idx = match[3]  # N9 in purine SMARTS pattern
            nucleobase_nitrogen_idx = n9_idx
            break
    elif has_pyrimidine:
        # In pyrimidines, N1 forms the glycosidic bond
        pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_pattern)
        for match in pyrimidine_matches:
            n1_idx = match[0]  # N1 in pyrimidine SMARTS pattern
            nucleobase_nitrogen_idx = n1_idx
            break

    if nucleobase_nitrogen_idx is None:
        return False, "Unable to identify nucleobase nitrogen for glycosidic bond"

    # Identify anomeric carbon of the sugar (connected to nucleobase)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    anomeric_carbon_idx = None
    for match in sugar_matches:
        # The anomeric carbon is the one connected to two oxygens
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
                if len(oxy_neighbors) == 2:
                    anomeric_carbon_idx = idx
                    break
        if anomeric_carbon_idx is not None:
            break

    if anomeric_carbon_idx is None:
        return False, "Anomeric carbon not found in sugar"

    # Check for N-glycosidic bond between nucleobase nitrogen and sugar anomeric carbon
    nucleobase_nitrogen = mol.GetAtomWithIdx(nucleobase_nitrogen_idx)
    anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)
    # Check if they are directly bonded
    is_bonded = mol.GetBondBetweenAtoms(nucleobase_nitrogen_idx, anomeric_carbon_idx)
    if not is_bonded:
        return False, "No N-glycosidic bond connecting nucleobase and sugar found"

    return True, "Molecule is a nucleoside with nucleobase attached to sugar via N-glycosidic bond"