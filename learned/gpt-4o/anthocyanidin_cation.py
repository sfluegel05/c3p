"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is characterized by a flavylium ion core with additional oxygenation.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Define a new SMARTS pattern for the flavylium core: [o+]1c(ccc2cc(O)cc(O)c12)c3ccccc3
    flavylium_pattern = Chem.MolFromSmarts("[o+]1cc(cc2cc(O)cc(O)c12)c3ccccc3")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return (False, "No flavylium core found or it is incorrectly structured")
    
    # Additional checks for correct oxygenation - There should be substantial substitution, typically hydroxyls
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    methoxy_pattern = Chem.MolFromSmarts("[OCH3]")
    ether_pattern = Chem.MolFromSmarts("[O]([#6])[#6]")  # Generic ether pattern

    # Count the oxygens bound directly to the core (hydroxyls/methoxy/etc.)
    oxy_group_matches = len(mol.GetSubstructMatches(hydroxyl_pattern)) + len(mol.GetSubstructMatches(methoxy_pattern)) + len(mol.GetSubstructMatches(ether_pattern))
    if oxy_group_matches < 2:
        return (False, "Insufficient oxygenation; at least two hydroxyl or ether groups expected")

    # Confirm presence of a positive charge on the oxygen atom
    positive_charge = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if not positive_charge:
        return (False, "Molecule lacks a positive charge indicating it is not a cation")

    return (True, "Contains a correctly oxygenated flavylium ion core with positive charge")