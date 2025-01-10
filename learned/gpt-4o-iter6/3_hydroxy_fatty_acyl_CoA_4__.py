"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4_(smiles: str):
    """
    Classifies a molecule as a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Pattern for 3-hydroxy group attached to a fatty acyl chain (-C-C(COC)O-)
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C;!R]")  # 3-hydroxy pattern
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy fatty acid moiety found"

    # Pattern for Coenzyme A via thioester linkage (-C(=O)SCCNC)
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)")  
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety via thioester linkage not found"

    # Pattern for negative charge in phosphate & diphosphate groups - ([O-])
    phosphates_pattern = Chem.MolFromSmarts("P([O-])([O-])")
    if len(mol.GetSubstructMatches(phosphates_pattern)) < 1:
        return False, "Deprotonated phosphate groups not found"

    return True, "The molecule is a 3-hydroxy fatty acyl-CoA(4-)"