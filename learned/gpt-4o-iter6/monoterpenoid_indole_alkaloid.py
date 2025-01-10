"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    Monoterpenoid indole alkaloids are complex structures derived from tryptophan and specialized 
    diisoprenoid precursors, with characteristic indole and expanded heterocyclic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced indole-like structures - adding hybrid and fused options
    indole_patterns = [
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3ccccc23"),  # Standard indole
        Chem.MolFromSmarts("C1=CC2=C3C=CC=CC3=NC2=C1"),  # Fused indoles
        Chem.MolFromSmarts("c1cc2cnc(c2cc1)N"),  # Indolinamine
        Chem.MolFromSmarts("c1cc2[nH]cnc2c1"),  # Pyrroloindole
        Chem.MolFromSmarts("n1c2ccc3c(c[nH]c3c2cc1)"),  # Complex indole patterns
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in indole_patterns):
        return False, "No recognizable indole-like structure"

    # Count nitrogen atoms with relaxation for tertiary amine or amide groups
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Insufficient nitrogen atoms for alkaloid complexity"

    # Simplified ring criteria to accommodate commonality in this class
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 4:
        return False, f"Not enough ring structures: found {num_rings}"

    # Check for the presence of common methoxy or ester groups (common in alkaloids)
    methoxy_pattern = Chem.MolFromSmarts("COC")  # Methoxy groups specifically
    ester_pattern = Chem.MolFromSmarts("C(=O)O")  # Useful ester moiety
    if not (mol.HasSubstructMatch(methoxy_pattern) or mol.HasSubstructMatch(ester_pattern)):
        return False, "Missing methoxy or ester groups"

    return True, "Contains indole core with necessary nitrogen and complex ring features typical of monoterpenoid indole alkaloids"