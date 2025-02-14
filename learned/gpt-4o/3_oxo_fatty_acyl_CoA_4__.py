"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class is characterized by the presence of a 3-oxo group and a Coenzyme A moiety.

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
    coa_smarts = "COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)N1C=NC2=C(N)N=CN=C12"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Check for 3-oxo group attached to a fatty acyl chain
    oxo_acyl_pattern = Chem.MolFromSmarts("C(=O)CC(=O)")
    if not mol.HasSubstructMatch(oxo_acyl_pattern):
        return False, "No 3-oxo group or incomplete acyl chain found"

    # Ensure there are long carbon chains indicating fatty acids
    chain_length = [atom.GetTotalNumHs() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    long_chain = any(cnum >= 14 for cnum in chain_length)
    if not long_chain:
        return False, "Carbon chain length too short for typical fatty acids"

    return True, "Valid 3-oxo-fatty acyl-CoA(4-) structure identified"

# Example usage:
smiles = "CCCCCCCCCCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles)
print(f"Result: {result}, Reason: {reason}")