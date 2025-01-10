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

    # Define a SMARTS pattern for monosaccharide rings (pyranose/furanose)
    # This pattern recognizes generic hexose rings (e.g., glucose-like, with oxygen in the ring)
    mono_pattern = Chem.MolFromSmarts('[C@@H]1([O])[C@@H]([C@@H]([C@@H]([C@H]([C@H]1O)O)O)O)O |$;@`,@`,@`,@`,@`,@`|')

    # Identify monosaccharide units in the molecule
    mono_matches = mol.GetSubstructMatches(mono_pattern)
    if len(mono_matches) < 2:
        return False, "Less than two monosaccharide units found"
    
    # Define a SMARTS pattern for the glycosidic bond (C-O-C linkage between two sugars using recursive SMARTS)
    # The pattern [C]-[O]-[C] is a simplistic way to identify glycosidic bonds
    glycosidic_bond_pattern = Chem.MolFromSmarts('[C]O[C]')

    # Identify glycosidic bonds in the molecule
    glyco_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glyco_matches) < 1:
        return False, "No glycosidic bonds found between sugar units"

    # Verify if the glycosidic bonds link two different monosaccharide units
    linked_monosaccharides = 0
    for bond_match in glyco_matches:
        # Check if the carbons in the glycosidic bond belong to different monosaccharide ring matches
        first_carbon_in_any_mono = any(bond_match[0] in mono for mono in mono_matches)
        second_carbon_in_any_mono = any(bond_match[2] in mono for mono in mono_matches)
        if first_carbon_in_any_mono and second_carbon_in_any_mono:
            linked_monosaccharides += 1

    # Verify exactly two sugar units are linked appropriately
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