"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class is characterized by the presence of a 3-oxo group and a Coenzyme A moiety,
    and a long fatty acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for Coenzyme A moiety
    # Coenzyme A pattern redesigned to be more flexible
    coa_pattern = Chem.MolFromSmarts("NC(=O)CC(COP(=O)([O-])OP(=O)([O-])OCC1O[C@H](n2cnc3c(N)ncnc23)[C@H]1O)O[C@H]1[C@H](O)[C@H](n2cnc3c(N)ncnc23)O[C@H]1COP(=O)([O-])OP([O-])(=O)OCC([C@@H](O)CNC(=O)CCNC(=O)C)[S]CC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Check for 3-oxo group attached to a fatty acyl chain
    # More precise pattern for 3-oxo group (O=C-C(=O))
    oxo_acyl_pattern = Chem.MolFromSmarts("C(=O)CC(=O)")
    if not mol.HasSubstructMatch(oxo_acyl_pattern):
        return False, "No 3-oxo group or incomplete acyl chain found"

    # Ensure there are long carbon chains indicating fatty acids
    carbon_chains = mol.GetSubstructMatches(Chem.MolFromSmarts("C" * 14))  # Look for at least 14 connected carbons
    if not carbon_chains:
        return False, "Carbon chain length too short for typical fatty acids"

    return True, "Valid 3-oxo-fatty acyl-CoA(4-) structure identified"

# Example usage:
smiles = "CCCCCCCCCCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles)
print(f"Result: {result}, Reason: {reason}")