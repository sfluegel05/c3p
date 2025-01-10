"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define purine and pyrimidine nucleobase patterns (more general)
    purine_base = Chem.MolFromSmarts('c1ncnc2ncnn12')  # purine ring system
    pyrimidine_base = Chem.MolFromSmarts('c1cncnc1')    # pyrimidine ring system

    # Check for nucleobase (purine or pyrimidine)
    has_nucleobase = mol.HasSubstructMatch(purine_base) or mol.HasSubstructMatch(pyrimidine_base)
    if not has_nucleobase:
        return False, "Nucleobase not found (purine or pyrimidine base)"

    # Define ribose and deoxyribose patterns (five-membered ring with oxygen)
    ribose = Chem.MolFromSmarts('C1[C@@H]([O])[C@H]([O])[C@@H](O)[C@H]1O')  # ribose sugar
    deoxyribose = Chem.MolFromSmarts('C1[C@@H]([O])[C@H]([O])[C@@H](O)[C@H]1')  # deoxyribose sugar

    # Check for sugar
    has_sugar = mol.HasSubstructMatch(ribose) or mol.HasSubstructMatch(deoxyribose)
    if not has_sugar:
        return False, "Sugar (ribose or deoxyribose) not found"

    # Define N-glycosidic bond between nucleobase and sugar (base attached to anomeric carbon via nitrogen)
    glycosidic_bond = Chem.MolFromSmarts('[n;!H0]-[C@H]1[C@@H]([O])[C@H]([O])[C@@H](O)[C@H]1[O]')  # N attached to anomeric carbon
    has_glycosidic_bond = mol.HasSubstructMatch(glycosidic_bond)
    if not has_glycosidic_bond:
        return False, "N-glycosidic bond between nucleobase and sugar not found"

    # Define phosphate group attached to 5' carbon of sugar
    phosphate = Chem.MolFromSmarts('OP(=O)(O)[O]')  # phosphate group
    phosphate_attachment = Chem.MolFromSmarts('[C@@H]1([O])[C@H]([O])[C@@H](O)[C@H](CO[P](=O)(O)[O])O1')  # sugar with phosphate at 5' position

    has_phosphate = mol.HasSubstructMatch(phosphate_attachment)
    if not has_phosphate:
        return False, "Phosphate group at 5' position not found"

    # Check for mono-, di-, tri-, or tetra-phosphate groups
    phosphate_counts = 0
    phosphates = Chem.MolFromSmarts('P(=O)(O)[O]')
    matches = mol.GetSubstructMatches(phosphates)
    phosphate_counts = len(matches)

    if phosphate_counts < 1 or phosphate_counts > 4:
        return False, f"Found {phosphate_counts} phosphate groups, but need between 1 and 4"

    return True, "Molecule is a nucleoside 5'-phosphate"

__metadata__ = {
    'chemical_class': {
        'name': "nucleoside 5'-phosphate",
        'definition': "A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated."
    }
}