"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for glycerophosphoserine backbone
    glycerophosphoserine_pattern = Chem.MolFromSmarts("OC[C@H](CON[P](=O)(O)OC[C@H](N)C(=O)O)O")
    if not mol.HasSubstructMatch(glycerophosphoserine_pattern):
        return False, "No glycerophosphoserine backbone found"

    # The acyl group pattern at the primary hydroxy position
    # Acyl group should be an ester link to the oxygen at the 1-position
    acyl_1_hydroxy_pattern = Chem.MolFromSmarts("OC(=O)C")
    # Check if the oxygen atom of glycerol is bonded to this pattern
    acyl_matches = mol.GetSubstructMatches(acyl_1_hydroxy_pattern)
    if not acyl_matches:
        return False, "No acyl group found at the 1-hydroxy position"
    
    for match in acyl_matches:
        # Verify that the acyl is correctly bonded in the context of the full backbone pattern
        acyl_bond = False
        for match1 in acyl_1_hydroxy_pattern.GetSubstructMatches(mol):
            if (mol.GetBondBetweenAtoms(match[0], match1[1]) or
                mol.GetBondBetweenAtoms(match[0], match1[0])):
                acyl_bond = True
                break
        
        if acyl_bond:
            return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"

    return False, "Acyl group not correctly found at 1-hydroxy position"