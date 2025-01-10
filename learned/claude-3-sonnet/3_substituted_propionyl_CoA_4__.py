"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA core structure components
    
    # Check for exactly 4 negative charges
    charge_pattern = Chem.MolFromSmarts("[O-]")
    charge_matches = mol.GetSubstructMatches(charge_pattern)
    if len(charge_matches) != 4:
        return False, f"Found {len(charge_matches)} negative charges, need exactly 4"

    # Check for adenine base
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No adenine base found"

    # Check for ribose with phosphate - more specific pattern
    ribose_phosphate = Chem.MolFromSmarts("OC1C(O)C(COP([O-])=O)OC(n2cnc3c(N)ncnc32)C1")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "No ribose-phosphate moiety found"

    # Check for pantetheine arm with thioester
    pantetheine_pattern = Chem.MolFromSmarts("NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "No pantetheine moiety found"

    # Check for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for diphosphate bridge
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)([O-])OP(=O)([O-])")
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "No diphosphate bridge found"

    # Verify acyl chain length (at least 3 carbons)
    acyl_chain = Chem.MolFromSmarts("SC(=O)CCC")
    if not mol.HasSubstructMatch(acyl_chain):
        return False, "Acyl chain too short or missing"

    # Count carbons and oxygens to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if p_count != 3:
        return False, "Must have exactly 3 phosphorus atoms"
    
    if c_count < 25:
        return False, "Insufficient carbon atoms for CoA structure"

    # All criteria met
    return True, "Contains CoA(4-) moiety with appropriate thioester linkage and substituted acyl chain"