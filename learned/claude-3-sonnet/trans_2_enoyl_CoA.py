"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: CHEBI:36415 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is an unsaturated fatty acyl-CoA resulting from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any 2,3-trans-enoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coenzyme A backbone pattern
    coa_pattern = Chem.MolFromSmarts("[N&R]1C2=NC(=NC2=NC1)N[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OC[C@@H]4[C@H]([C@@H]([C@@H](O4)N5C=NC6=C5N=CN=C6N)O)OP(=O)(O)O)O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found"

    # Look for trans-2-enoyl group pattern
    enoyl_pattern = Chem.MolFromSmarts("[CX3](/C=C/[CX3](=O))")
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No trans-2-enoyl group found"

    # Look for fatty acid chain (minimum 6 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Check molecular weight - trans-2-enoyl-CoA typically 800-1500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800 or mol_wt > 1500:
        return False, "Molecular weight outside typical range for trans-2-enoyl-CoA"

    return True, "Contains coenzyme A backbone with trans-2-enoyl fatty acid chain"