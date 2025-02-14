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

    # Look for glycerol backbone pattern (C(O)C(O)C with allowances for stereochemistry)
    # Here each carbon atom is bonded to an oxygen atom, typical of triglycerides
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for 3 ester groups (-C(=O)O-)
    # We want exactly three ester groups for a triglyceride
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, expected 3"

    # Check for chains to ensure they are part of fatty acids
    # We define a minimal pattern for fatty acid chains to allow flexibility
    long_chain_min_length = Chem.MolFromSmarts("CCCC")
    fatty_acid_matches = mol.GetSubstructMatches(long_chain_min_length)
    
    # Ensure total chains longer than three carbons are found at least three times
    if len(fatty_acid_matches) < 3:
        return False, f"Missing sufficient fatty acid-like chains, detected {len(fatty_acid_matches)}"

    # Molecular weight check
    # Formerly used as a hard check, now more of a guidepost. For diverse triglycerides, MW can vary but typically above 350 Da
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 350:
        return False, "Molecular weight is quite low, usually indicative of triglyceride variants or shorter chains"

    return True, "SMILES represents a triglyceride with a glycerol backbone and three ester groups"