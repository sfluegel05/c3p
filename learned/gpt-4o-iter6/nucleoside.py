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
    
    # Expanded nucleobase pattern - consider common forms and tautomers
    nucleobase_patterns = [
        Chem.MolFromSmarts("c1ncnc2ncnc12"),  # Purine, common form
        Chem.MolFromSmarts("c1[nH]cnc2ncnc12"),  # Purine variant
        Chem.MolFromSmarts("c1nc([nH])c(=O)[nH]c1=O"),  # Uracil/Thymine variant
        Chem.MolFromSmarts("c1[nH]c(=O)[nH]c(=O)[nH]1"),  # Uracil/Thymine
        Chem.MolFromSmarts("c1nc2[nH]c[nH]c2[nH]1"),  # Adenine
        Chem.MolFromSmarts("c1nc2c([nH]1)ncnc2=O"),  # Guanine
        Chem.MolFromSmarts("c1[nH]c(=O)ncn1"),  # Cytosine
        Chem.MolFromSmarts("n1cc[nH]c1=O"),  # Any pyrimidinone variant, generic
        Chem.MolFromSmarts("n1cnc2ncnc2c1=O"),  # Hypoxanthine
        Chem.MolFromSmarts("n1cnc2ncnc2c1"),  # Xanthine
    ]
    
    # Check for nucleobase substructures
    if not any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns):
        return False, "No nucleobase found"
    
    # Ribose and deoxyribose sugar patterns with flexible stereochemistry
    ribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@H]([C@@H]1O)O)C")
    deoxyribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H]([C@H]1)O)C")
    sugar_patterns = [ribose_pattern, deoxyribose_pattern]
    
    # Check for sugar substructures
    if not any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns):
        return False, "No ribose or deoxyribose sugar found"
    
    # Linkage check between sugar and nucleobase
    nucleoside_linkage = Chem.MolFromSmarts("[n,c]1:n2-[c,nH]=,:[nH,c.o]-c3[o,n]c([o,n]c([o,n]c3-o2-c4c([o,n]c([o,n]c4-o)))o)n1")
    if not mol.HasSubstructMatch(nucleoside_linkage):
        return False, "Missing correct linkage between nucleobase and sugar"
    
    # If all checks pass, it's a nucleoside
    return True, "Contains a nucleobase and a ribose or deoxyribose sugar with appropriate linkage"