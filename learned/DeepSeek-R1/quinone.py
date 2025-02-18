"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone has a fully conjugated cyclic dione structure derived from aromatic compounds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all ketone groups (C=O)
    ketone_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[#6]=O'))
    if len(ketone_matches) < 2:
        return False, f"Found {len(ketone_matches)} ketone groups, need at least 2"
    
    # Check for conjugated cyclic dione structure using SMARTS pattern
    # Matches two ketones in a conjugated 6-membered ring system
    quinone_pattern = Chem.MolFromSmarts('[O]=C1C=CC(=O)C=C1')
    if mol.HasSubstructMatch(quinone_pattern):
        return True, "Contains conjugated cyclic dione structure (basic quinone)"
    
    # More general pattern for conjugated cyclic diones (any ring size)
    conjugated_dione = Chem.MolFromSmarts('[O]=C1C=CC(=O)C=C1')
    if mol.HasSubstructMatch(conjugated_dione):
        return True, "Contains conjugated cyclic dione structure"