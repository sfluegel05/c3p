"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a sulfonate group (-SO3-) attached to a carbon.
    The carbon at position 1 can be attached to hydrogens, a carbon chain, or other groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for basic alkanesulfonate pattern:
    # [#6] - any carbon
    # -[SX4] - connected to sulfur with 4 connections
    # (=[OX1])(=[OX1]) - two double bonded oxygens
    # [O-] - one negatively charged oxygen
    alkanesulfonate_pattern = Chem.MolFromSmarts(
        '[#6]-[SX4](=[OX1])(=[OX1])[O-]'
    )
    
    if not mol.HasSubstructMatch(alkanesulfonate_pattern):
        return False, "No alkanesulfonate group found"

    # Get matches
    matches = mol.GetSubstructMatches(alkanesulfonate_pattern)
    
    # Check each match
    for match in matches:
        carbon_idx = match[0]
        sulfur_idx = match[1]
        
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Verify carbon is bonded to sulfur with single bond
        bond = mol.GetBondBetweenAtoms(carbon_idx, sulfur_idx)
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        return True, "Contains alkanesulfonate group (C-SO3-)"

    return False, "No valid alkanesulfonate group found"