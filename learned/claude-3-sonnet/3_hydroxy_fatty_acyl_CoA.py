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

    # Check for CoA core structure
    # Adenine + ribose + phosphates pattern
    coa_pattern = Chem.MolFromSmarts("[nX2r5:1]1c[nX2r5]c2c1nc([nX2r5]c2N)N")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found (missing adenine)"

    # Check for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for 3-hydroxy pattern (-CH2-CH(OH)-CH2-)
    # This pattern should be adjacent to the thioester
    hydroxy_pattern = Chem.MolFromSmarts("[CX4H2]-[CX4H]([OX2H])-[CX4H2]-C(=O)[SX2]")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy group found adjacent to thioester"

    # Check for pantetheine arm pattern
    pantetheine_pattern = Chem.MolFromSmarts("CC(C)(COP(O)OP)[CH](O)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing or incomplete pantetheine arm"

    # Check for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    phosphate_matches = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_matches < 3:
        return False, f"Insufficient phosphate groups (found {phosphate_matches}, need at least 3)"

    # Count carbons in the fatty acid chain
    # First, find the thioester carbon
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

    return True, "Contains CoA moiety, thioester linkage, 3-hydroxy group, and fatty acid chain"