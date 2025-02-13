"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is characterized by a glycerol backbone with two acyl chains and a phosphocholine head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with sn configuration [C@H] or [C@@H]
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with correct stereochemistry not found"

    # Look for phosphocholine head group
    phosphocholine_pattern = Chem.MolFromSmarts("P(=O)(OCC[N+](C)(C)C)")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine head group not found"
        
    # Check for two acyl chains (esters connected to glycerol backbone)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H](COC(=O))")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} acyl chains, need exactly 2"

    # Verify the stereochemistry if further needed - (e.g., 'sn' configuration)
    stereochemistry_correct = mol.HasSubstructMatch(Chem.MolFromSmarts("O[C@H](CO)"))
    if not stereochemistry_correct:
        return False, "Incorrect stereochemistry, not an sn-glycerol derivative"
    
    return True, "Contains glycerol backbone with two acyl chains and a phosphocholine head group"