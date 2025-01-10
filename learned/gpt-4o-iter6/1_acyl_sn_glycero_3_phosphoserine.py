"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    This classification is based on the presence of glycerophosphoserine backbone and an acyl
    substituent at the 1-hydroxy position, attached through an ester linkage.

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

    # Define a more generic SMARTS pattern for the glycerophosphoserine backbone
    glycerophosphoserine_pattern = Chem.MolFromSmarts("P(OC[C@H](O)CO)(OC[C@H](N)C(=O)O)O")
    
    # Alter to allow variable regions and account for common stereo centers in this class
    if not mol.HasSubstructMatch(glycerophosphoserine_pattern):
        return False, "No glycerophosphoserine backbone found"

    # Define a SMARTS pattern for the acyl group linked through an ester at the primary position
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)O[C@H](O)C")
    acyl_alternative_pattern = Chem.MolFromSmarts("C(=O)OCC(O)C")  # Allows for absence of explicit stereochemistry
    if not mol.HasSubstructMatch(acyl_chain_pattern) and not mol.HasSubstructMatch(acyl_alternative_pattern):
        return False, "No acyl group found specifically at the 1-hydroxy position"

    # If both patterns are found, this is a positive match
    return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"