"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: CHEBI:15892 aliphatic alcohol

An aliphatic alcohol is an alcohol derived from an aliphatic compound.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alcohol group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1,OX2H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if not alcohol_matches:
        return False, "No alcohol group found"

    # Check for aliphatic carbon chain attached to the alcohol group
    aliphatic_pattern = Chem.MolFromSmarts("[CX4H3,CX4H2]~[CX4H2]~[CX4H2]")
    for match in alcohol_matches:
        alcohol_atom = mol.GetAtomWithIdx(match)
        if any(neighbor.HasSubstructMatch(aliphatic_pattern) for neighbor in alcohol_atom.GetNeighbors()):
            break
    else:
        return False, "No aliphatic carbon chain attached to alcohol group"

    # Count rotatable bonds to check for long aliphatic chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Insufficient rotatable bonds for aliphatic chains"

    return True, "Contains an alcohol group attached to an aliphatic carbon chain"