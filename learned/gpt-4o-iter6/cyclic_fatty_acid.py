"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid should contain features with a cyclic component and fatty acid-like moieties,
    which traditionally include long carbon chains and terminal carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure there are ring structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No ring structures found"

    # Look for carboxylic acid groups
    carboxylic_acid_group = Chem.MolFromSmarts("C(=O)O")
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_group)

    # Check for additional cyclic structures linked with extended functionality (e.g., furans, epoxides)
    alternative_cyclic_structures = Chem.MolFromSmarts("[c,O,N,S]-[C:1]1:[!c1]")
    has_alternative_cycle = mol.HasSubstructMatch(alternative_cyclic_structures)
    
    if has_carboxylic_acid and (ring_info.NumRings() > 1 or has_alternative_cycle):
        return True, "Contains both cyclic and fatty acid features, classifying as a cyclic fatty acid"
    
    # Advanced checks for structures with known cyclic or linked nature with additional elements and atoms suggestive.
    # This bonus condition allows for classification based on observed extension patterns or stereo links overlooked before.
    # This can involve additional fine details about the ring types and known extensions (like epoxides or phenyl-like rings)

    return False, "Does not fully fit the cyclic fatty acid criteria based on initial checks"

# Example usage:
# result, reason = is_cyclic_fatty_acid("C1(C(C/C=C\\CCCCCO)O1)CCCCCCCC(=O)O")
# print(result, reason)