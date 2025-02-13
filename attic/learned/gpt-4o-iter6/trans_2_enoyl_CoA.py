"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # SMARTS pattern for trans double bond between second and third carbon in fatty acid chain
    trans_double_bond_pattern = Chem.MolFromSmarts("C\\C=C\\C")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No trans double bond found between second and third carbon"

    # SMARTS pattern for thioester linkage (C(=O)SCC)
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage not found"

    # Check for the presence of key structural components indicative of coenzyme A
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)OC[C@@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not detected"

    return True, "Molecule matches all patterns for trans-2-enoyl-CoA"