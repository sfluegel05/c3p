"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime contains an aliphatic chain and an aldehyde group
    converted to an oxime (R-C=NOH), without presence in aromatic systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for the classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for oxime group pattern "C=N-O"
    oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX1H1]")
    oxime_matches = mol.GetSubstructMatches(oxime_pattern)
    if not oxime_matches:
        return False, "No oxime group found"

    # Ensure oxime is part of an aliphatic environment
    for match in oxime_matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        nitrogen_atom = mol.GetAtomWithIdx(match[1])

        # Oxime group should not be part of a ring or aromatic system
        if carbon_atom.IsInRing() or nitrogen_atom.IsInRing():
            return False, "Oxime group is part of a ring"

        # Check neighboring atoms to ensure they are not part of an aromatic system
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return False, "Oxime group is part of an aromatic system"
    
        # If we reached here, the oxime is in an aliphatic setting
        return True, "Contains an aliphatic aldoxime group"

    return False, "Oxime group is not part of a purely aliphatic chain"