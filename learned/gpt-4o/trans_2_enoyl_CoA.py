"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    
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

    # Check for coenzyme A adenine component
    coa_pattern = Chem.MolFromSmarts("n1cnc2c(n1)c[nH]c2")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Not a CoA molecule (missing adenine base)"

    # Check for thiol ester linkage to CoA: -C(=O)SCCNC(=O)
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Missing thiol ester linkage"

    # Check for trans double bond at the 2,3 position: C\C=C\
    trans_double_bond_pattern = Chem.MolFromSmarts("C\\C=C\\")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No 2,3-trans double bond found"
    
    return True, "Molecule fits the trans-2-enoyl-CoA structure criteria"