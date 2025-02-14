"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is a lipid containing one or more hydroperoxy (-OOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for hydroperoxy group (-OOH or -OO-)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX1]")  # Matches -OOH and -OO-
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy group found"

    # 2. Check for sufficient heavy atoms to be a lipid
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_heavy_atoms < 10:
         return False, f"Too few heavy atoms ({num_heavy_atoms}), not likely a lipid"

    # 3. Check for at least one carbon atom (to handle edge cases)
    carbon_pattern = Chem.MolFromSmarts("[CX4,CX3]")
    if not mol.HasSubstructMatch(carbon_pattern):
        return False, "Not a carbon based molecule"

    return True, "Contains at least one hydroperoxy group and is likely a lipid"