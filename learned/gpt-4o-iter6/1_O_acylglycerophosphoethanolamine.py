"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
from rdkit import Chem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    A 1-O-acylglycerophosphoethanolamine has a glycerophosphoethanolamine structure with an O-acyl substituent at the 1-position of the glycerol fragment.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern checking for 1-O-acylglycerol moiety with stereochemistry
    o_acyl_pattern = Chem.MolFromSmarts("[C@@H](COC(=O)[C])[O][C@H](O)CO")
    if not mol.HasSubstructMatch(o_acyl_pattern):
        return False, "No O-acyl ester correctly positioned at the 1-position with necessary stereochemistry"

    # Improved pattern for phosphoethanolamine group, targeting precise connectivity and format
    phosphoethanolamine_pattern = Chem.MolFromSmarts("P(=O)([O-])OC[C@H](N)C")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No correctly oriented phosphoethanolamine group present"

    # Check for appropriate number of glycerol ester linkages and phosphate presence
    ester_linkages_pattern = Chem.MolFromSmarts("COC(=O)C")
    phosphate_count = sum(1 for _ in mol.GetSubstructMatches(ester_linkages_pattern))
    if phosphate_count < 1:
        return False, f"Insufficient glycerophosphoester linkages, found {phosphate_count}"
    
    return True, "Contains glycerophosphoethanolamine structure with O-acyl substituent at the 1-position"