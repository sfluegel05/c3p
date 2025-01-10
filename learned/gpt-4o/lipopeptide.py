"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide consists of a peptide with attached lipid moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define SMARTS pattern for a peptide bond (-CO-NH-)
    peptide_bond_pattern = Chem.MolFromSmarts("[C](=O)[N]")
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"

    # Enhanced lipid moieties detection
    # Extended to consider long chains, branching, and cyclic hydrocarbons
    lipid_patterns = [
        Chem.MolFromSmarts("[CH2,CH3]{8,}"),  # Simple long carbon chains with 8+ atoms
        Chem.MolFromSmarts("C1CCCCCCCCCC1"),  # Simple cyclic structures
        Chem.MolFromSmarts("C~C(~C)~C"),         # Branching patterns
    ]

    lipid_like_found = any(mol.HasSubstructMatch(pattern) for pattern in lipid_patterns)
    if not lipid_like_found:
        return False, "No lipid-like structures found"

    # Check for cyclic peptides
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings > 0:
        return True, "Contains cyclic structures, peptide bonds, and lipid moieties"

    return True, "Contains peptide bonds and lipid-like moieties"

# Example usage:
# smiles_string = "CCCCC(=O)NCC(=O)NCCC(=O)OC1=C(NC2=C1C=CC=C2)O"
# is_lipopeptide(smiles_string)