"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
from rdkit import Chem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Classifies a molecule as trans-2-enoyl-CoA based on its SMILES string.
    Checks for a trans double bond, thioester linkage, and Coenzyme A moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as trans-2-enoyl-CoA
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for trans double bond connected to a thioester linkage
    trans_double_bond_pattern = Chem.MolFromSmarts("C/C=C\\C(=O)S")
    if not mol.HasSubstructMatch(trans_double_bond_pattern):
        return False, "No suitable trans double bond with thioester found"
    
    # Individually define key features of Coenzyme A and search for them collectively
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)[C@H](O)C(C)(C)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine component not detected"
    
    phosphoadenosine_pattern = Chem.MolFromSmarts("OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(=O)O)n2cnc3c(N)ncnc23")
    if not mol.HasSubstructMatch(phosphoadenosine_pattern):
        return False, "Phosphoadenosine component not detected"
    
    # Combine components to verify Coenzyme A structure as a whole
    if not (mol.HasSubstructMatch(pantetheine_pattern) and mol.HasSubstructMatch(phosphoadenosine_pattern)):
        return False, "Coenzyme A component not fully detected"
    
    return True, "Matches key features of trans-2-enoyl-CoA"