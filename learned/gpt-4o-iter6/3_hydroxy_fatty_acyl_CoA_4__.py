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

    # Pattern for 3-hydroxy group attached to a chiral center in fatty acyl chain
    hydroxy_chiral_pattern = Chem.MolFromSmarts("[C@@H](O)[CH2]")
    if not mol.HasSubstructMatch(hydroxy_chiral_pattern):
        return False, "No chiral 3-hydroxy group found on fatty acid chain"

    # Pattern for Coenzyme A thioester linkage with a specific context
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)C[C@H](O)C(C)(C)COP([O-])([O-])=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety with thioester linkage not found"

    # Pattern for at least one deprotonated phosphate group in CoA
    phosphate_pattern = Chem.MolFromSmarts("P([O-])([O-])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "At least two deprotonated phosphate groups not found"

    # Confirm presence of fatty acid chain length and unsaturation handling (exemplarily)
    chain_length_pattern = Chem.MolFromSmarts("CCCC(CCCCCC)CCCC")
    if not mol.HasSubstructMatch(chain_length_pattern):
        return False, "Appropriate long fatty acyl chain not found"

    return True, "The molecule is a 3-hydroxy fatty acyl-CoA(4-)"