"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is a glycerol backbone with three fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone, allow for stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for 3 ester groups (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, expected 3"

    # Check for long carbon chains indicative of fatty acids
    # A simple approach is to count successive carbon connections, a more thorough computational check can be employed
    fatty_acid_chain_pattern = Chem.MolFromSmarts("CCCC")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_chain_pattern)
    if len(fatty_acid_matches) < 3:
        return False, f"Missing long fatty acid chains, got {len(fatty_acid_matches)}"

    # Molecular weight check remains but less strict
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 500:
        return False, "Molecular weight too low for traditional triglyceride but could be a shorter variant"

    return True, "SMILES represents a triglyceride with a glycerol backbone and three ester groups"