"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: acyl-CoA compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed between coenzyme A and a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester linkage (-C(=O)S-) - this is essential for acyl-CoA
    thioester_patterns = [
        Chem.MolFromSmarts("[CX3](=[OX1])[SX2]"),
        Chem.MolFromSmarts("C(=O)S")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in thioester_patterns):
        return False, "No thioester linkage found"

    # Check for adenine base with flexible connection point
    adenine_patterns = [
        Chem.MolFromSmarts("c1nc(N)c2ncnc2n1"),
        Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in adenine_patterns):
        return False, "No adenine moiety found"

    # Check for phosphate groups with flexible patterns
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)O"),
        Chem.MolFromSmarts("OP([O,OH])(=O)[O,OH]"),
        Chem.MolFromSmarts("OP(O)(O)=O")
    ]
    phosphate_count = 0
    for pattern in phosphate_patterns:
        if pattern:
            phosphate_count += len(mol.GetSubstructMatches(pattern))
    if phosphate_count < 2:  # Allow for some variation in phosphate representation
        return False, f"Insufficient phosphate groups found"

    # Check for pantetheine portion with flexible pattern
    pantetheine_patterns = [
        Chem.MolFromSmarts("NCCC(=O)NCCS"),
        Chem.MolFromSmarts("NCCSC(=O)"),
        Chem.MolFromSmarts("NC(=O)CCNC(=O)")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in pantetheine_patterns):
        return False, "No pantetheine moiety found"

    # Check for ribose sugar with very flexible pattern
    ribose_patterns = [
        Chem.MolFromSmarts("OC1C(O)C(O)C(O)C1O"),
        Chem.MolFromSmarts("OC1CCCO1"),  # Basic furanose pattern
        Chem.MolFromSmarts("C1OC(CO)C(O)C1"),
        Chem.MolFromSmarts("C1OC(COP)C(O)C1")  # Connected to phosphate
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in ribose_patterns):
        return False, "No ribose sugar found"

    # Check for characteristic CoA features with flexible patterns
    coa_patterns = [
        Chem.MolFromSmarts("CC(C)(COP)"),  # Geminal dimethyl
        Chem.MolFromSmarts("SCCNC(=O)"),   # Core connection
        Chem.MolFromSmarts("COP(O)OP")     # Phosphate linkage
    ]
    if not all(mol.HasSubstructMatch(pattern) for pattern in coa_patterns):
        return False, "Missing key CoA structural features"

    # If all essential elements are present, it's an acyl-CoA
    return True, "Contains complete CoA moiety with thioester linkage to acyl group"