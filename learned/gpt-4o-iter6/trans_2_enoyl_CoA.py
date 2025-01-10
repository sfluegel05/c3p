"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    An enoyl-CoA compound characterized by the formal condensation of a thioester linkage with a trans 2,3-enoic acid.

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

    # Generalized SMARTS for trans 2,3-enoic acid motif.
    trans_double_bond_pattern = Chem.MolFromSmarts("[C;!R][C;!R]=[C;!R][C;!R]")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No trans double bond found in a plausible position"

    # SMARTS for thioester linkage (O=C-SC-) allowing flexible chain lengths
    thioester_pattern = Chem.MolFromSmarts("*C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Relaxed SMARTS for key CoA moiety structure.
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)O[C@H]1O[C@H]([C@@H](O)[C@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A component not detected"

    if mol.GetNumAtoms() > 40:
        return False, "Unexpectedly complex; verify SMILES conforms to a plausible trans-2-enoyl-CoA"

    return True, "Molecule matches all critical patterns for trans-2-enoyl-CoA"