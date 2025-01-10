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

    # Check for generic O-acyl ester group at the 1-position of glycerol fragment, assuming no specific stereochemistry
    o_acyl_pattern = Chem.MolFromSmarts("[C](OC(=O))[C](O)CO")
    if not mol.HasSubstructMatch(o_acyl_pattern):
        return False, "No acyl ester at the 1-position of glycerol fragment"

    # Check for phosphoethanolamine group (using a simplified, more generic pattern P(=O)(O)OCCN)
    phosphoethanolamine_pattern = Chem.MolFromSmarts("P(=O)(O)OCCN")
    if not mol.HasSubstructMatch(phosphoethanolamine_pattern):
        return False, "No phosphoethanolamine group attached"

    return True, "Contains glycerophosphoethanolamine structure with O-acyl substituent at the 1-position"