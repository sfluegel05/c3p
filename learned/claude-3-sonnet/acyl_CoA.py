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

    # Check for adenine base (part of CoA)
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine moiety found (required for CoA)"

    # Check for phosphate groups (CoA has 3 phosphates)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, f"Found only {len(phosphate_matches)} phosphate groups, need at least 3"

    # Check for thioester linkage (-C(=O)S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for pantetheine portion (part of CoA)
    # Looking for NCCC(=O)NCCS pattern
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine moiety found"

    # Check for ribose sugar (part of CoA)
    ribose_pattern = Chem.MolFromSmarts("OC1C(O)C(O)C(O)C1")
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No ribose sugar found"

    # Additional check for characteristic geminal dimethyl group in CoA
    geminal_dimethyl = Chem.MolFromSmarts("CC(C)(COP)")
    if not mol.HasSubstructMatch(geminal_dimethyl):
        return False, "Missing characteristic geminal dimethyl group"

    # If all structural elements are present, it's likely an acyl-CoA
    return True, "Contains CoA moiety with thioester linkage to acyl group"