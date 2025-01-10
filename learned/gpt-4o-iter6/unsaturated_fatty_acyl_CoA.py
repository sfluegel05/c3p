"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA is defined as a fatty acyl-CoA with a coenzyme A moiety linked
    via a thioester bond to an unsaturated fatty acid (contains C=C bonds).

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

    # Define functional group patterns:
    # Look for coenzyme A moiety pattern (specific parts)
    coenzyme_a_pattern = Chem.MolFromSmarts("P(=O)(O)OC[C@H]1O[C@H](CN2C=NC3=C2N=CN=C3N)[C@H](O)[C@H]1OP(O)(=O)O")
    if not mol.HasSubstructMatch(coenzyme_a_pattern):
        return False, "Coenzyme A moiety not found"

    # Look for thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for carbon-carbon double bond (unsaturation)
    c_double_bond_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(c_double_bond_pattern):
        return False, "No carbon-carbon double bond found, not unsaturated"

    # Check for sufficient carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Carbon chain too short, found only {carbon_count} carbons"

    return True, "Molecule is an unsaturated fatty acyl-CoA with required functional groups"

# Example usage:
smiles_example = "CCC\C=C\C=C\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)N1C=NC2=C1N=CN=C2N"
is_unsaturated_fatty_acyl_CoA(smiles_example)