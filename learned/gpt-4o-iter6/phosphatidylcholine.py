"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerophosphocholine with two acyl substituents.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the glycerophosphocholine backbone
    phosphocholine_pattern = Chem.MolFromSmarts("OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No glycerophosphocholine backbone found"
    
    # Look for two acyl groups (ester linkages attached to glycerol)
    acyl_pattern = Chem.MolFromSmarts("OC(=O)C")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl groups, need at least 2"
    
    # Check for appropriate attachments and stereochemistry at glycerol
    glycerol_structure = Chem.MolFromSmarts("[C@H](O[*])COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(glycerol_structure):
        return False, "Glycerol center or sterochemistry not matched"

    # If all checks pass, classify as phosphatidylcholine
    return True, "Contains glycerophosphocholine backbone with two acyl substituents"