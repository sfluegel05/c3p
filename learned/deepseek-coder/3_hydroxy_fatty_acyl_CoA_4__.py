"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA(4-) is a molecule with a 3-hydroxy fatty acid chain attached to a CoA moiety,
    with deprotonated phosphate and diphosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for 3-hydroxy fatty acid chain
    # Pattern: at least 4 carbons with a hydroxyl group at position 3
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4]([OH])[CX4]C(=O)")
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No 3-hydroxy fatty acid chain found"

    # Check charge state - should have 4 negative charges
    negative_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    if negative_charges != 4:
        return False, f"Expected 4 negative charges, found {negative_charges}"

    # Check molecular weight - should be >700 Da for typical 3-hydroxy fatty acyl-CoA
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for 3-hydroxy fatty acyl-CoA"

    # Count carbons in fatty acid chain - should be at least 4
    fatty_acid_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if fatty_acid_carbons < 14:  # CoA + minimum fatty acid chain
        return False, f"Too few carbons ({fatty_acid_carbons}) for 3-hydroxy fatty acyl-CoA"

    return True, "Contains CoA moiety with 3-hydroxy fatty acid chain and 4 negative charges"