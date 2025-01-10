"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA has a coenzyme A moiety linked via a thioester to an
    unsaturated fatty acid chain (with one or more C=C bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coenzyme A moiety pattern (ribose, phosphates, and adenine)
    ribose_phosphate_adenine_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H]([C@H](O)[C@@H]1)[OP](=[O])([O])[O][C@H][n1]c2[n][cH][n][c]12")
    if not mol.HasSubstructMatch(ribose_phosphate_adenine_pattern):
        return False, "No coenzyme A moiety found"

    # Look for thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"
        
    # Look for unsaturation - examine presence of carbon-carbon double bond(s)
    c_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(c_double_bond_pattern):
        return False, "No C=C double bond found, not unsaturated"

    # Check for a long enough carbon chain (typical fatty acids have 12+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Carbon chain too short, only {carbon_count} carbons"
    
    return True, "Molecule contains CoA, thioester linkage, C=C double bond, and sufficient carbon chain"

# Example usage:
smiles_example = "CCC\C=C\C=C\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"
is_unsaturated_fatty_acyl_CoA(smiles_example)