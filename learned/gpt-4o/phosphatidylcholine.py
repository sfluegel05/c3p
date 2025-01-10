"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
from rdkit import Chem

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

    # Glycerol backbone pattern with correct stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[C@H](OC(=O)[CX4,CX3])[C@@H](OC(=O)[CX4,CX3])COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone and/or stereochemistry not found"

    # Phosphocholine head group pattern without specific charges
    phosphocholine_pattern = Chem.MolFromSmarts("OP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "Phosphocholine head group not found"

    # Check for two ester-linked acyl chains connected to glycerol
    ester_pattern = Chem.MolFromSmarts("OC(=O)[CX4,CX3]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester-linked acyl chains, need exactly 2"

    return True, "Contains glycerol backbone with two ester-linked acyl chains and a phosphocholine head group"