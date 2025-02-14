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
    
    # Check for aldehyde or alcohol group at one end
    end_group_pattern = Chem.MolFromSmarts("[CH2O,CH2C=O]")
    has_end_group = mol.HasSubstructMatch(end_group_pattern)
    
    # Check for pentose backbone (5 carbons with 4 oxygens attached)
    pentose_pattern = Chem.MolFromSmarts("[CR2][OR2][CR2][OR2][CR2][OR2][CR2][CR2]")
    has_pentose_backbone = mol.HasSubstructMatch(pentose_pattern)
    
    # Classify as aldopentose if it has an aldehyde/alcohol group at one end and a pentose backbone
    if has_end_group and has_pentose_backbone:
        return True, "Contains a pentose backbone with an aldehyde or alcohol group at one end"
    else:
        return False, "Does not meet the criteria for an aldopentose"