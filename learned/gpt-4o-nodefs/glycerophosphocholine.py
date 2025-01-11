"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine has a glycerol backbone typically connected via ester or ether linkages
    and a characteristic phosphocholine group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with ester or ether linkages
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COP(*)=O)(OC*)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ester or ether linkages found"
        
    # Look for phosphocholine group
    # The pattern should include a phosphate group attached to a quaternary ammonium ion
    phosphocholine_pattern = Chem.MolFromSmarts("[P](=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    return True, "Contains glycerol backbone with ester/ether linkages and a phosphocholine group"