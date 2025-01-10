"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl at position n (omega).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the longest carbon chain
    lwc = Chem.rdmolops.GetLongestPath(mol)
    chain_length = len(lwc)
    if chain_length < 6:
        return False, f"Carbon chain too short has only {chain_length} atoms"

    # Check if the first atom in this chain is part of a carboxyl group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches or lwc[0] not in next(zip(*carboxyl_matches)):
        return False, "No carboxyl group at position 1 of the main chain"

    # Check if the last atom in this chain is part of a terminal hydroxyl group
    if mol.GetAtomWithIdx(lwc[-1]).GetAtomicNum() != 8:
        return False, "No omega-hydroxyl group at end of the main chain"

    return True, "Contains carboxyl group at position 1 and omega-hydroxyl group at end of the main chain as per omega-hydroxy fatty acid definition."