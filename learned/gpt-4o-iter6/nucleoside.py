"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    
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
    
    # Expanded nucleobase patterns to cover major nucleobases
    nucleobase_patterns = [
        Chem.MolFromSmarts("c1ncnc2ncnc12"),           # Purine-like base
        Chem.MolFromSmarts("n1cnc2[nH]cnc2[nH]c1"),    # Adenine
        Chem.MolFromSmarts("c1nc2c([nH]1)ncnc2=O"),    # Guanine
        Chem.MolFromSmarts("c1c[nH]cnc1=O"),           # Thymine and Uracil
        Chem.MolFromSmarts("n1ccnc2[nH]cnc12")         # Cytosine-like 
    ]
    
    # Check for nucleobase substructures
    contains_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    
    if not contains_nucleobase:
        return False, "No nucleobase found"
    
    # Ribose and deoxyribose sugar patterns with flexible stereochemistry
    ribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@H]([C@@H]1O)O)C")        # Ribose pattern
    deoxyribose_pattern = Chem.MolFromSmarts("CC1COC(OCC1)O")                        # Deoxyribose pattern
    
    # Check for sugar substructures
    contains_sugar = mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not contains_sugar:
        return False, "No ribose or deoxyribose sugar found"
    
    # Generalize N-glycosidic linkage pattern
    n_glycosidic_pattern = Chem.MolFromSmarts("[OD2]C1OCC(O1)[ND]")  # N-glycosidic linkage
    has_glycosidic_bond = mol.HasSubstructMatch(n_glycosidic_pattern)
    
    if not has_glycosidic_bond:
        return False, "Missing N-glycosidic linkage between sugar and nucleobase"
    
    return True, "Contains a nucleobase and a ribose or deoxyribose sugar with appropriate linkage"