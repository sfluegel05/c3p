"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine has a glycerol backbone potentially with ester or ether-linked fatty acids
    and a phosphocholine group.

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

    # Look for glycerol backbone with potential ester/ether bonds
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](CO*)(O*)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with ether or ester linkages found"
        
    # Look for phosphocholine group
    # Relaxed pattern to consider different orientations/conformations
    phosphocholine_pattern = Chem.MolFromSmarts("[P](=O)(O)OC[C@@H]([N+](C)(C)C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    return True, "Contains both glycerol backbone with fatty acids or ethers and a phosphocholine group"