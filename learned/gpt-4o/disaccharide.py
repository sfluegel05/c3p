"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound with two monosaccharides joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Common pattern for detecting hexose and pentose rings
    pyran_furan_pattern = Chem.MolFromSmarts('[C@@H]1O[C@H](CO)[C@@H](O)[C@@H](O)[C@H]1O')
    
    # Identify monosaccharide (pyranose or furanose) units in the molecule
    mono_matches = mol.GetSubstructMatches(pyran_furan_pattern)
    if len(mono_matches) < 2:
        return False, "Less than two monosaccharide-like rings found"

    # General pattern for glycosidic bonds (C-O-C linkage)
    glycosidic_bond_pattern = Chem.MolFromSmarts('COC')

    # Identify glycosidic bonds in the molecule
    glyco_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glyco_matches) < 1:
        return False, "No glycosidic bonds found between monosaccharide rings"

    # Verify if the glycosidic bonds link two distinct monosaccharide units
    linked_monosaccharides = 0
    for bond_match in glyco_matches:
        # Iterate through glycosidic bond found (COC identified)
        for mono in mono_matches: 
            if (bond_match[0] in mono) and (bond_match[2] in mono):
                linked_monosaccharides += 1
                break

    # There needs to be at least one glycosidic bond linking two rings
    if linked_monosaccharides >= 1:
        return True, "Contains two monosaccharide units joined by a glycosidic bond"
    else:
        return False, "Insufficient linkage between monosaccharide units"

# Test examples
examples = [
    "O([C@H]1[C@H](O)[C@H](O)C(O[C@@H]1CO)O)[C@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO",
    "O([C@H]1[C@@H](O)[C@H](O[C@@H](O)[C@@H]1O)CO)[C@@H]2OC[C@H](O)[C@H](O)[C@H]2O"
]

for example in examples:
    result, reason = is_disaccharide(example)
    print(f"SMILES: {example} -> Is disaccharide? {result}. Reason: {reason}")