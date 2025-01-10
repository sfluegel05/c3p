"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    The classification is based on having a glycerophosphoserine backbone with an acyl group 
    at the 1-hydroxy position, linked through an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern to capture 1-acyl-sn-glycero-3-phosphoserine structure, focusing more on acyl link
    glycerophosphoserine_pattern = Chem.MolFromSmarts("P(OC[C@H](O)CO)(OC[C@H](N)C(=O)O)O")
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)OC[C@H](O)CO")
    
    # Check for specific backbone
    if not mol.HasSubstructMatch(glycerophosphoserine_pattern):
        return False, "No glycerophosphoserine backbone found"

    # Check for specific acyl linkage
    if not mol.HasSubstructMatch(acyl_chain_pattern):
        return False, "No acyl group at 1-hydroxy position found"

    return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"