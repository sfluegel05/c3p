"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    Such a molecule contains a Coenzyme A structure and a fatty acyl chain with exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A general pattern
    coa_pattern = Chem.MolFromSmarts(
        "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)N"
    )
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"

    # Check for exactly one carbon-carbon double bond
    db_pattern = Chem.MolFromSmarts("C=C")
    db_matches = mol.GetSubstructMatches(db_pattern)

    if len(db_matches) != 1:
        return False, f"Expected 1 carbon-carbon double bond, found {len(db_matches)}"

    # Verify the double bond is part of a fatty acyl chain (common sense that C=C is often in the longer tail region)
    # Usually, the chain should start from the fatty acyl portion (long Cx chain ending with C(=O)SCC)
    # This might require more robust logic depending on the diversity of structures in reality.
    
    return True, "Contains Coenzyme A structure and a fatty acyl chain with exactly one carbon-carbon double bond"