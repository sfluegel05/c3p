"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is characterized by a glycerol backbone with two acyl chains
    and a phosphocholine head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone with flexible stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with correct stereochemistry not found"

    # Phosphocholine head group pattern
    phosphocholine_pattern = Chem.MolFromSmarts("O=P([O-])(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine head group not found"

    # Check for two acyl chains connected to glycerol. Assuming any length fatty acid chains.
    acyl_chain_pattern = Chem.MolFromSmarts("OC(=O)C")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) < 2:
        return False, f"Found {len(acyl_matches)} acyl chains, need at least 2"

    return True, "Contains glycerol backbone with two acyl chains and a phosphocholine head group"