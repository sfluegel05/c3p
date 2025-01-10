"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organometalloidal_compound(smiles: str):
    """
    Determine if a molecule is an organometalloidal compound based on its SMILES string.
    Organometalloidal compounds have bonds between metalloids (e.g., As) and organic carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for different patterns indicative of organometalloidal compounds
    # Direct As-C bonds, part of an organic group
    arsenic_carbon_pattern = Chem.MolFromSmarts("[As]-[C]")
    # Organyl groups with functional arsenic
    organyl_arsenic_pattern = Chem.MolFromSmarts("[C][As](=O)(O)")
    # Expanded pattern with aromatic or complex bonds
    aromatic_arsenic_pattern = Chem.MolFromSmarts("c[As]c")

    # Checking all defined patterns
    if mol.HasSubstructMatch(arsenic_carbon_pattern):
        return True, "Contains arsenic-carbon bonds indicative of organometalloidal compounds"
    elif mol.HasSubstructMatch(organyl_arsenic_pattern):
        return True, "Contains functional arsenic organyl group"
    elif mol.HasSubstructMatch(aromatic_arsenic_pattern):
        return True, "Contains aromatic complex with arsenic characteristic of organometalloidal compounds"

    return False, "Does not meet organometalloidal structural criteria"