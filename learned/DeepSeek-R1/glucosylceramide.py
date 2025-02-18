"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: CHEBI glucosylceramide (cerebrosides with glucose head group)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    Glucosylceramides consist of a ceramide (sphingosine + fatty acid) linked to a glucose molecule via beta-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for glucose (beta-D-glucopyranose) substructure
    glucose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@H](O[C@@H]1CO)O)O)O)O")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Check glycosidic bond (glucose connected via oxygen to ceramide)
    glycosidic_o = Chem.MolFromSmarts("[O;X2][C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(glycosidic_o):
        return False, "No glycosidic oxygen bond found"

    # Check for ceramide part: sphingosine (long chain with NH and OH) + fatty acid amide
    # Look for amide group connected to a long aliphatic chain
    amide_pattern = Chem.MolFromSmarts("[NH1X3][C]=O")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"

    # Check that the amide is part of a long chain (fatty acid)
    # Assuming fatty acid has at least 12 carbons (common in ceramides)
    fatty_acid_length = 12
    for match in amide_matches:
        amide_n = mol.GetAtomWithIdx(match[0])
        neighbor = amide_n.GetNeighbors()[0]
        chain = Chem.MolFromSmarts(f"[{neighbor.GetSymbol()}]-[#6](=O)-[#6](-[#6]){'{'}{fatty_acid_length-2},}")
        if mol.HasSubstructMatch(chain):
            break
    else:
        return False, "Fatty acid chain too short or not found"

    # Check sphingosine part: long chain with hydroxyl and amino groups
    # Simplified pattern: look for a chain with at least one hydroxyl and the amide-connected NH
    sphingosine_pattern = Chem.MolFromSmarts("[C]([OH])[C]([NH])")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "Sphingosine structure not found"

    return True, "Contains beta-D-glucose linked via glycosidic bond to ceramide"