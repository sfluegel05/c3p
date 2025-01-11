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

    # Define nucleobase patterns
    adenine = Chem.MolFromSmarts('n1cnc2c1ncnc2')  # adenine
    guanine = Chem.MolFromSmarts('n1[cH]nc2c1ncnc2O')  # guanine
    cytosine = Chem.MolFromSmarts('n1ccnc1=O')  # cytosine
    thymine = Chem.MolFromSmarts('Cc1c[nH]c(=O)[nH]c1=O')  # thymine
    uracil = Chem.MolFromSmarts('O=C1NC=CC(=O)N1')  # uracil

    nucleobases = [adenine, guanine, cytosine, thymine, uracil]

    # Define ribose and deoxyribose sugar patterns (five-membered ring with oxygen)
    ribose = Chem.MolFromSmarts('[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](CO)O1')  # ribose
    deoxyribose = Chem.MolFromSmarts('C1(O)C(O)C(O)C(CO)O1')  # deoxyribose

    sugars = [ribose, deoxyribose]

    # Define phosphate group patterns (mono-, di-, tri-, tetra-phosphate)
    monophosphate = Chem.MolFromSmarts('OP(=O)(O)O')  # monophosphate
    diphosphate = Chem.MolFromSmarts('OP(=O)(O)OP(=O)(O)O')  # diphosphate
    triphosphate = Chem.MolFromSmarts('OP(=O)(O)OP(=O)(O)OP(=O)(O)O')  # triphosphate
    tetraphosphate = Chem.MolFromSmarts('OP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)O')  # tetraphosphate

    phosphates = [monophosphate, diphosphate, triphosphate, tetraphosphate]

    # Check for nucleobase attached to sugar via N-glycosidic bond
    nucleoside_found = False
    for base in nucleobases:
        for sugar in sugars:
            # Combine base and sugar patterns with an anomeric bond
            nucleoside_pattern = Chem.MolFromSmarts('[nH]1c2c(nc[nH]2)ncn1[C@H]3O[C@H](CO)[C@@H](O)[C@H]3O')
            if mol.HasSubstructMatch(nucleoside_pattern):
                nucleoside_found = True
                break
        if nucleoside_found:
            break

    if not nucleoside_found:
        return False, "Nucleobase attached to sugar not found"

    # Check for sugar with phosphate group at 5' position
    phosphate_found = False
    for phosphate in phosphates:
        # Define sugar-phosphate pattern: phosphate attached to 5' carbon of sugar
        sugar_phosphate_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](OP(=O)(O)O)CO1')  # ribose-5'-phosphate
        if mol.HasSubstructMatch(sugar_phosphate_pattern):
            phosphate_found = True
            break
        # Check deoxyribose as well
        deoxy_sugar_phosphate_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](OP(=O)(O)O)CO1')  # deoxyribose-5'-phosphate
        if mol.HasSubstructMatch(deoxy_sugar_phosphate_pattern):
            phosphate_found = True
            break

    if not phosphate_found:
        return False, "Phosphate group at 5' position not found"

    return True, "Molecule is a nucleoside 5'-phosphate"

__metadata__ = {
    'chemical_class': {
        'name': "nucleoside 5'-phosphate",
        'definition': "A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated."
    }
}