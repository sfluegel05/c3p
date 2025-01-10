"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is defined as a polyhydroxy aldehyde or its intramolecular hemiacetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for aldehyde group (-C(=O)H)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No terminal aldehyde group found"

    # Check for multiple hydroxyl groups (polyhydroxy structure)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient number of hydroxyl groups for polyhydroxy structure"

    # Check for potential cyclic form of sugar (hemiacetal)
    ring_info = mol.GetRingInfo()
    ring_shapes = ring_info.AtomRings()
    if not any(len(ring) in [5, 6] for ring in ring_shapes):
        return False, "No 5 or 6-membered rings detected for cyclic form"

    # Verify each identified ring has Ether (O-C-O) linkage and at least additional 2 hydroxyl groups
    ether_and_oh_count = 0
    for ring in ring_shapes:
        if len(ring) in [5, 6]:
            ethers = [atom for atom in ring if mol.GetAtomWithIdx(atom).GetSymbol() == 'O']
            if len(ethers) == 1:
                ether_and_oh_count += sum(1 for atom in ethers if mol.GetAtomWithIdx(atom).GetDegree() == 2)
    
    if ether_and_oh_count < 1:
        return False, "Cyclic form doesn't meet ether linkage and hydroxyl condition for an aldose"

    return True, "Contains aldehyde group with polyhydroxy and potential cyclic structure matching an aldose"