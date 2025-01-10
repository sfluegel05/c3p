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
        return None, "Invalid SMILES string"

    # Improved Coenzyme A general pattern
    coa_pattern = Chem.MolFromSmarts(
        "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)C"  # Shortened for explanation, customize as needed
        "OP(=O)(O)O"
        "P(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)N"
    )
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"
    
    # Check for exactly one carbon-carbon double bond
    db_pattern = Chem.MolFromSmarts("C=C")
    db_matches = mol.GetSubstructMatches(db_pattern)

    if len(db_matches) != 1:
        return False, f"Expected 1 carbon-carbon double bond, found {len(db_matches)}"

    # Ensure the double bond is part of the fatty acyl chain
    for match in db_matches:
        atom_indices = {atom.GetIdx() for atom in mol.GetAtoms()}
        if all(idx in atom_indices for idx in match):
            return True, "Contains Coenzyme A structure and a fatty acyl chain with exactly one carbon-carbon double bond"
    
    return False, "Double bond not part of a valid fatty acyl chain"