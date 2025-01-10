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
    
    # Extended nucleobase patterns
    purine_patterns = [
        Chem.MolFromSmarts("ncnc2ncnc2"), # Basic purine scaffold
        Chem.MolFromSmarts("ncnc2ncn[nH]2"), # Including tautomers and modifications
        Chem.MolFromSmarts("n1[nH]cnc2ncnc12") # Ring closure variations
    ]
    
    pyrimidine_patterns = [
        Chem.MolFromSmarts("c1cnc[nH]c1=O"), # Basic pyrimidine scaffold
        Chem.MolFromSmarts("c1c[nH]cnc1=O"), # Including tautomers and modifications
        Chem.MolFromSmarts("c1c[nH]c(=O)n(c1=O)") # Oxy groups variations
    ]
    
    # Check for purine or pyrimidine nucleobase-like substructures
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in purine_patterns + pyrimidine_patterns)
    
    if not has_nucleobase:
        return False, "No nucleobase found"
    
    # Sugar pattern development
    ribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](CO)C1") # ribose with stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](C)C1") # deoxyribose
    
    # Check for ribose or deoxyribose-like substructures
    has_sugar = mol.HasSubstructMatch(ribose_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"
    
    # Verify connectivity between bases and sugars: looking for N-glycosidic bond
    # Well-connected sugar and nucleobase is the final verification
    n_glycosidic_bond = Chem.MolFromSmarts("OC1C(O)C(O)COC1n") # Check for bond to nitrogen on nucleobase
    has_glycosidic_bond = mol.HasSubstructMatch(n_glycosidic_bond)
    
    if not has_glycosidic_bond:
        return False, "Missing N-glycosidic linkage between sugar and nucleobase"
    
    return True, "Contains a nucleobase and a ribose or deoxyribose sugar with appropriate linkage"

# Example SMILES testing, should test against a list of known nucleoside SMILES for validation