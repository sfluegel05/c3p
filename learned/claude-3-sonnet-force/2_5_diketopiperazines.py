"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:34614 2,5-diketopiperazine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_5_diketopiperazine(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is a cyclic compound containing a piperazine-2,5-dione skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for piperazine-2,5-dione core as part of a cyclic structure
    core_pattern = Chem.MolFromSmarts("C1NC(=O)CN(C)C1=O")
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "No piperazine-2,5-dione core found"

    # Check for appropriate substituents on the piperazine ring
    substituent_pattern = Chem.MolFromSmarts("[C;R]1[N;R](C(=O)N[C;R]([C;R](=O)[C;R]([N;R]1[C;R]))[C;R])([C;R])[C;R]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    if not substituent_matches:
        return False, "Substituents on piperazine ring are not typical of 2,5-diketopiperazines"

    # Consider stereochemistry, if applicable
    try:
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception:
        pass  # Ignore stereochemistry if conformational analysis fails

    # Additional checks for specific structural patterns or substructures
    # ...

    return True, "Contains a cyclic piperazine-2,5-dione core with appropriate substituents"