"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA core structure - multiple patterns to catch different representations
    # Adenine patterns
    adenine_patterns = [
        Chem.MolFromSmarts("n1cnc2c(N)ncnc12"),  # Common form
        Chem.MolFromSmarts("N1C=NC2=C1N=CN=C2N"),  # Alternative form
        Chem.MolFromSmarts("[nX2r5]1c[nX2r5]c2c1[nX2r5]c([nX2r5]c2N)N")  # More specific form
    ]
    
    has_adenine = any(mol.HasSubstructMatch(pattern) for pattern in adenine_patterns)
    if not has_adenine:
        return False, "No adenine moiety found in CoA structure"

    # Check for ribose-phosphate part
    ribose_phosphate = Chem.MolFromSmarts("OCC1OC(n2cnc3c2ncnc3N)C(O)C1OP(O)(O)=O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing ribose-phosphate structure"

    # Check for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for 3-hydroxy pattern with possible stereochemistry
    hydroxy_patterns = [
        # R configuration
        Chem.MolFromSmarts("[CX4H2]-[CX4H]([OX2H])-[CX4H2]-C(=O)[SX2]"),
        # S configuration
        Chem.MolFromSmarts("[CX4H2]-[CX4@H]([OX2H])-[CX4H2]-C(=O)[SX2]"),
        Chem.MolFromSmarts("[CX4H2]-[CX4@@H]([OX2H])-[CX4H2]-C(=O)[SX2]")
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in hydroxy_patterns):
        return False, "No 3-hydroxy group found adjacent to thioester"

    # Check for pantetheine arm
    pantetheine_patterns = [
        Chem.MolFromSmarts("CC(C)(COP(O)OP)C(O)C(=O)NCCC(=O)NCCS"),
        Chem.MolFromSmarts("CC(C)(COP([O-])OP)C(O)C(=O)NCCC(=O)NCCS")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in pantetheine_patterns):
        return False, "Missing or incomplete pantetheine arm"

    # Check for phosphate groups
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(O)(=O)O"),
        Chem.MolFromSmarts("OP([O-])(=O)O")
    ]
    
    total_phosphates = sum(len(mol.GetSubstructMatches(pattern)) 
                          for pattern in phosphate_patterns)
    if total_phosphates < 3:
        return False, f"Insufficient phosphate groups (found {total_phosphates}, need at least 3)"

    # Count carbons in the fatty acid chain
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if thioester_matches:
        thioester_carbon = thioester_matches[0][0]
        # Count carbons in the fatty acid portion
        fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]")
        fatty_acid_carbons = len([atom for atom in mol.GetAtoms() 
                                if atom.GetAtomicNum() == 6 and 
                                atom.GetDegree() <= 4 and
                                len(Chem.GetShortestPath(mol, atom.GetIdx(), thioester_carbon)) <= 20])
        
        if fatty_acid_carbons < 4:
            return False, f"Fatty acid chain too short (found {fatty_acid_carbons} carbons)"

    return True, "Contains CoA moiety with adenine, ribose-phosphate, pantetheine arm, thioester linkage, 3-hydroxy group, and appropriate fatty acid chain"