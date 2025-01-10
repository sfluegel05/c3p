"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    It checks for a trans double bond with a thioester attached to the CoA moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for trans double bond (E configuration) near start of chain
    trans_double_bond_pattern = Chem.MolFromSmarts("C/C=C\\C(=O)")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No trans double bond found in plausible position"

    # SMARTS for thioester linkage (O=C-SC-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Expanded SMARTS pattern for generic Coenzyme A moiety recognizing essential parts
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A component not detected"

    # Check for reasonable complexity (atom count)
    if mol.GetNumAtoms() > 150:  
        return False, "Unexpectedly complex structure; verify that this is a likely trans-2-enoyl-CoA"

    return True, "Molecule matches all critical patterns for trans-2-enoyl-CoA"