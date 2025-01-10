"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23553 diol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound containing exactly two hydroxyl (-OH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First find all potential OH groups
    # [OX2H1] matches oxygen with 2 connections, 1 hydrogen
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Early exit if less than 2 OH groups
    if len(oh_matches) < 2:
        return False, f"Found only {len(oh_matches)} hydroxyl groups, need 2"

    # Find carboxylic acids and exclude their OH groups
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    carboxyl_ohs = {match[2] for match in carboxyl_matches}  # Get OH oxygen indices

    # Find enols and exclude their OH groups
    enol_pattern = Chem.MolFromSmarts("[CX3]=[CX3][OX2H1]")
    enol_matches = mol.GetSubstructMatches(enol_pattern)
    enol_ohs = {match[2] for match in enol_matches}  # Get OH oxygen indices

    # Find phenols and exclude their OH groups
    phenol_pattern = Chem.MolFromSmarts("c[OX2H1]")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    phenol_ohs = {match[1] for match in phenol_matches}  # Get OH oxygen indices

    # Get valid alcoholic OH groups (exclude carboxylic acids, enols, and phenols)
    valid_oh_matches = []
    for match in oh_matches:
        oh_idx = match[0]
        if oh_idx not in carboxyl_ohs and oh_idx not in enol_ohs and oh_idx not in phenol_ohs:
            oh_oxygen = mol.GetAtomWithIdx(oh_idx)
            # Check that OH is connected to carbon
            for neighbor in oh_oxygen.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    if not neighbor.GetIsAromatic():  # Not aromatic carbon
                        valid_oh_matches.append(match)
                        break

    # Check if we have exactly 2 valid OH groups
    if len(valid_oh_matches) != 2:
        return False, f"Found {len(valid_oh_matches)} valid alcoholic hydroxyl groups, need exactly 2"

    # Get the carbons that the OH groups are attached to
    oh_carbons = []
    for match in valid_oh_matches:
        oh_oxygen = mol.GetAtomWithIdx(match[0])
        for neighbor in oh_oxygen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                oh_carbons.append(neighbor)

    # Verify both OH groups are attached to sp3 carbons
    if not all(carbon.GetHybridization() == Chem.HybridizationType.SP3 for carbon in oh_carbons):
        return False, "Not all hydroxyl groups are attached to sp3 carbons"

    return True, "Contains exactly two alcoholic hydroxyl groups"