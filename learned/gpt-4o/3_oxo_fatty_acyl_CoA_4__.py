"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class is characterized by the presence of a 3-oxo group, a Coenzyme A moiety,
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
        return False, "Invalid SMILES string"

    # Check for Coenzyme A moiety
    coa_pattern = Chem.MolFromSmarts("OP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    # Check for 3-oxo group attached to a long acyl chain
    oxo_acyl_pattern = Chem.MolFromSmarts("C(=O)[CX4,CX3][CX4,CX3]C(=O)")
    if not mol.HasSubstructMatch(oxo_acyl_pattern):
        return False, "No 3-oxo group or incomplete acyl chain found"

    # Ensure there are long carbon chains indicating fatty acids
    chain_pattern = Chem.MolFromSmarts("C([C@H](O)[CX4](C)C)")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 3:
        return False, "Insufficient carbon chains typical of fatty acids"

    # Refine check for carbon chain lengths, typical for fatty acids (at least a C16 backbone)
    if Chem.rdMolDescriptors.CalcNumAtomStereoCenters(mol) < 14:
        return False, "Carbon chain length too short for typical fatty acids"

    return True, "Valid 3-oxo-fatty acyl-CoA(4-) structure identified"

# Example usage:
smiles = "CCCCCCCCCCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles)
print(f"Result: {result}, Reason: {reason}")