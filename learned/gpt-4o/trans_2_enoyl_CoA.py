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

    # Use a broader SMARTS pattern for the adenine component of CoA
    coa_adenine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)N")  # Improved adenine pattern
    if not mol.HasSubstructMatch(coa_adenine_pattern):
        return False, "Not a CoA molecule, adenine not found"

    # Identify the linkage of coenzyme A from a fatty acyl chain through a thioester linkage
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)SC")  # Simplified to allow variable length
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Missing thioester linkage"

    # Specifically identify the trans double bond at the second position, more specific to trans geometry
    trans_double_bond_pattern = Chem.MolFromSmarts("C/C=C/C")
    double_bond_matches = mol.GetSubstructMatches(trans_double_bond_pattern)
    if len(double_bond_matches) == 0:
        return False, "No 2,3-trans double bond found"

    # Ensure all key CoA features are present
    coa_full_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C1O[C@@H](COP(O)=O)C[C@H](O1)n2cnc3c(ncnc23)N")
    if not mol.HasSubstructMatch(coa_full_pattern):
        return False, "Full CoA structure not confirmed"

    return True, "Molecule fits the trans-2-enoyl-CoA structure criteria"