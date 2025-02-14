"""
Classifies: CHEBI:18379 nitrile
"""
from rdkit import Chem

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile is characterized by a carbon atom triply bonded to a nitrogen atom (C#N),
    where the carbon atom is directly substituted to a non-hydrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a nitrile group.
    nitrile_pattern = Chem.MolFromSmarts("[CX4;!H0]#[N]")

    # Find matches of the nitrile pattern
    matches = mol.GetSubstructMatches(nitrile_pattern)
    if not matches:
        return False, "Does not contain a nitrile group (C#N) where the carbon is not attached to a hydrogen atom"

    # Check for at least one non-hydrogen atom attached to the carbon. We can check the neighbor atoms
    for match in matches:
        nitrile_carbon_index = match[0]
        nitrile_carbon = mol.GetAtomWithIdx(nitrile_carbon_index)
        
        #If all the neighbors are hydrogens, then this is not a proper nitrile.
        is_valid_nitrile = False
        for neighbor in nitrile_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() != 1:
                is_valid_nitrile = True
                break
        if not is_valid_nitrile:
           return False, "Does not contain a nitrile group (C#N) where the carbon is not attached to a non-hydrogen atom"

    return True, "Contains a nitrile group (C#N) where the carbon is not attached to hydrogen(s)."