"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: CHEBI:57926 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    A trans-2-enoyl-CoA is an unsaturated fatty acyl-CoA that results from the formal condensation
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

    # Look for CoA backbone pattern (more flexible pattern)
    coa_pattern = Chem.MolFromSmarts("[C;R](C)(C)(CO[P@@](=O)(O)O[P@@](=O)(O)OC[C@H]1[C@@H]([C@@H]([C@H](O1)n2cnc3c2ncnc3N)O)OP(=O)(O)O)C(=O)NCCCC(=O)NCCSC(=O)"
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # Look for trans double bond at position 2 (more flexible pattern)
    trans_pattern = Chem.MolFromSmarts("C\C=C/C(=O)[CX4]")
    if not mol.HasSubstructMatch(trans_pattern):
        return False, "No trans double bond at position 2 found"

    # Look for fatty acid chain (long carbon chain attached to the carbonyl)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "No fatty acid chain found"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Fatty acid chain too short"

    # Check molecular weight - CoA derivatives typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for CoA derivative"

    return True, "Contains CoA backbone and trans double bond at position 2, with fatty acid chain attached"