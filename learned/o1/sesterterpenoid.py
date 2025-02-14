"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from sesterterpenes composed of five isoprene units (C25 skeleton),
    possibly modified by rearrangement or loss of small groups (e.g., methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 30:
        return False, f"Carbon count is {c_count}, not typical for sesterterpenoids (approx. 25 carbons)"
    
    # Check for isoprene units
    # Isoprene unit: C=C-C-C=C
    # We will look for multiple occurrences of this pattern
    isoprene_pattern = Chem.MolFromSmarts('C(=C)CC=C')
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 3:
        return False, f"Found {len(isoprene_matches)} isoprene units, less than expected for sesterterpenoids"
    
    # Check for oxygen-containing functional groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, unlikely to be a terpenoid (which are oxygenated terpenes)"
    
    # Check for typical terpenoid functional groups
    # e.g., hydroxyl, carbonyl, ester, ether
    functional_groups = [
        Chem.MolFromSmarts('[OX2H]'),   # Hydroxyl group
        Chem.MolFromSmarts('C=O'),      # Carbonyl group
        Chem.MolFromSmarts('COC'),      # Ether group
        Chem.MolFromSmarts('OC=O'),     # Ester group
    ]
    fg_found = False
    for fg in functional_groups:
        if mol.HasSubstructMatch(fg):
            fg_found = True
            break
    if not fg_found:
        return False, "No typical terpenoid functional groups found"

    return True, "Molecule meets the criteria for a sesterterpenoid"