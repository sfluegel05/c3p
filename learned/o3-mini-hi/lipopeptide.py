"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: Lipopeptide - a compound consisting of a peptide with attached lipid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide consists of a peptide component (evident by amide bonds) 
    and a lipid component (manifested as a long aliphatic chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is recognized as a lipopeptide, False otherwise.
        str: A reason detailing the classification.
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide (amide) bonds
    # Typical peptide/amide bond fragment: -N-C(=O)-
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    peptide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(peptide_matches) < 1:
        return False, "No amide (peptide) bonds found; no peptide component detected"

    # Check for lipid component: we look for a long uninterrupted aliphatic chain.
    # Here we use a simple heuristic pattern "CCCCCCCC" to match at least 8 contiguous sp3 carbons.
    lipid_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No long aliphatic chain found; lipid component missing"

    # Optional additional checks (e.g., molecular weight, chain length) could be performed here.
    # For now, if both substructures are detected, we classify this as a lipopeptide.
    return True, "Contains both a peptide component (amide bonds) and a lipid component (long aliphatic chain)"