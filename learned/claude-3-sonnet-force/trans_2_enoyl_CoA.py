"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: CHEBI:36356 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS
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
    
    # Check for trans-2-enoyl pattern
    trans_enoyl_pattern = Chem.MolFromSmarts("[CX3]([H])=C(/[CH3,CH2,CH1,C])\\C(=O)")
    if not mol.HasSubstructMatch(trans_enoyl_pattern):
        return False, "No trans-2-enoyl pattern found"
    
    # Check for CoA backbone
    coa_pattern = rdFMCS.FindMoleculeChemicalEnvironment(mol, Chem.MolFromSmiles("C(C(C(=O)NCCC(=O)NCCS))O)OP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc23)C(O)C1OP(=O)(O)O"))
    if coa_pattern is None:
        return False, "No CoA backbone found"
    
    # Check for long carbon chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Carbon chain too short"
    
    # Additional checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    
    if c_count < 20:
        return False, "Too few carbons for trans-2-enoyl-CoA"
    if o_count < 10:
        return False, "Too few oxygens for trans-2-enoyl-CoA"
    if n_count != 9:
        return False, "Incorrect number of nitrogens for trans-2-enoyl-CoA"
    if p_count != 3:
        return False, "Incorrect number of phosphorus atoms for trans-2-enoyl-CoA"
    
    return True, "Contains trans-2-enoyl pattern and CoA backbone"