"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is an N-oxide where there are three organic groups bonded to the nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for N-oxide using SMARTS. Check that the nitrogen is bonded to 3 carbons and an oxygen.
    # Also check that each carbon is connected to another carbon.
    n_oxide_pattern = Chem.MolFromSmarts("[N+](-[O-])([CX4])([CX4])[CX4]")
    if not mol.HasSubstructMatch(n_oxide_pattern):
        return False, "Not a tertiary amine N-oxide"
    
    #Check if carbons are indeed part of an organic group by checking they have at least one other neighbour that is not H or a heteroatom.
    match = mol.GetSubstructMatch(n_oxide_pattern)
    nitrogen_atom = mol.GetAtomWithIdx(match[0])
    carbon_atoms = [mol.GetAtomWithIdx(match[i]) for i in range(2,5)]
    for carbon_atom in carbon_atoms:
        neighbors = carbon_atom.GetNeighbors()
        organic_neighbor = False
        for neighbor in neighbors:
             if neighbor.GetIdx() != nitrogen_atom.GetIdx(): # Avoid checking back towards the central N atom
                if neighbor.GetAtomicNum() not in [1, 7, 8, 16, 15]: # Check for H and common heteroatoms
                    organic_neighbor = True
        if not organic_neighbor:
             return False, "Nitrogen is not bonded to three organic groups"

    return True, "Tertiary amine oxide detected"