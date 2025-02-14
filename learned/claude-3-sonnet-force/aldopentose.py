"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:18237 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a (potential) aldehyde group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for aldehyde group (both [CH]=O and C=O forms)
    aldehyde_pattern = Chem.MolFromSmarts("[CH]=O,C=O")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)
    
    # Check for pentose backbone (5 carbons with 4 oxygens attached)
    pentose_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CHX4][CHX4][CHX3]")
    has_pentose_backbone = mol.HasSubstructMatch(pentose_pattern)
    
    # Check for 4 hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    has_4_hydroxyls = len(hydroxyl_matches) == 4
    
    # Classify as aldopentose if it has an aldehyde group, a pentose backbone, and 4 hydroxyl groups
    if has_aldehyde and has_pentose_backbone and has_4_hydroxyls:
        return True, "Contains a pentose backbone with an aldehyde group and 4 hydroxyl groups"
    else:
        return False, "Does not meet the criteria for an aldopentose"