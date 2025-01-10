"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a given SMILES string corresponds to a polyprenol phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect isoprene units; allow for various stereochemistries using more general unsaturation patterns
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C=C")  # Simplified and generalized isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 3:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 3 for a polyprenol chain"

    # Detect phosphate group; capture both mono- and diphosphate groups, protonated and deprotonated
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")  # Cover more forms including fully protonated
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate or diphosphate group found"

    # Check for the linkage (ester bond) between polyprenol and phosphate group
    linkage_pattern = Chem.MolFromSmarts("C-O-P(=O)(O)O")  # Check for C-O-P linkage
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "Polyprenol chain not appropriately connected to phosphate group"

    return True, "Molecule is a polyprenol phosphate"

# Example usage:
# smiles = "CC(C)=CCC\C(C)=C\CC\C(C)=C\CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/CC\C(C)=C/COP(O)(=O)OP(O)(O)=O"
# print(is_polyprenol_phosphate(smiles))