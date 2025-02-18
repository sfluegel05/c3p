"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride consists of a glycerol backbone with three fatty acid chains attached via ester bonds.

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

    # Check for exactly three ester groups (-O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 3"

    # Check glycerol backbone with three ester attachments
    glycerol_pattern = Chem.MolFromSmarts(
        "[CH2]([OX2][CX3](=[OX1]))[CH]([OX2][CX3](=[OX1]))[CH2]([OX2][CX3](=[OX1]))"
    )
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Missing glycerol backbone with three ester attachments"

    # Verify fatty acid chain length via molecular weight (typical triglycerides > 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:  # Lowered threshold to accommodate some natural variations
        return False, f"Molecular weight too low ({mol_wt:.1f} < 400)"

    # Count total carbons (typical triglycerides have >45 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 30:
        return False, f"Insufficient carbons ({carbon_count} < 30)"

    return True, "Glycerol backbone with three fatty acid esters meeting criteria"