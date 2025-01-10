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

    # Check for glycerol backbone with O-acyl at 1-position ([C@@H](COC(=O))
    glycerol_acyl_pattern = Chem.MolFromSmarts("[C@@H](COC(=O))[C@@H](O)CO")
    acyl_match = mol.HasSubstructMatch(glycerol_acyl_pattern)
    if not acyl_match:
        return False, "No acyl ester at the 1-position of glycerol fragment"

    # Check for phosphoethanolamine group (P(O)(=O)OCCN)
    phosphoethanolamine_pattern = Chem.MolFromSmarts("POC[C@@H](O)COC(=O)CP(=O)(O)OCCN")
    phospho_match = mol.HasSubstructMatch(phosphoethanolamine_pattern)
    if not phospho_match:
        return False, "No phosphoethanolamine group attached"

    return True, "Contains glycerophosphoethanolamine structure with O-acyl substituent at the 1-position"