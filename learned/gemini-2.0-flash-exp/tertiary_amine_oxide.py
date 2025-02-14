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
    n_oxide_pattern = Chem.MolFromSmarts("[N+](-[O-])([C])([C])[C]")
    if not mol.HasSubstructMatch(n_oxide_pattern):
          return False, "Not an N-oxide with three carbons."

    #Check if the nitrogen is positively charged.
    match = mol.GetSubstructMatch(n_oxide_pattern)
    nitrogen_atom = mol.GetAtomWithIdx(match[0])
    if nitrogen_atom.GetFormalCharge() != 1:
         return False, "Nitrogen is not positively charged"
    
    # Check that the nitrogen is connected only to the oxygen and three carbons, that are not hydrogens
    n_neighbors = nitrogen_atom.GetNeighbors()
    if len(n_neighbors) != 4:
          return False, "Nitrogen is not connected to 3 carbons and 1 oxygen"
    oxygen_count = 0
    carbon_count = 0
    for neighbor in n_neighbors:
        if neighbor.GetAtomicNum() == 8:
            oxygen_count += 1
        elif neighbor.GetAtomicNum() == 6:
            carbon_count += 1
        elif neighbor.GetAtomicNum() == 1:
             return False, "Nitrogen is connected to hydrogens"
    if oxygen_count != 1 or carbon_count != 3:
         return False, "Nitrogen is not connected to 3 carbons and 1 oxygen"
        
    # Check if carbons are indeed part of an organic group by checking they have at least one other neighbour that is not H or a heteroatom.
    
    carbon_atoms = [mol.GetAtomWithIdx(match[i]) for i in range(2,5)]
    for carbon_atom in carbon_atoms:
        neighbors = carbon_atom.GetNeighbors()
        organic_neighbor = False
        for neighbor in neighbors:
             if neighbor.GetIdx() != nitrogen_atom.GetIdx(): # Avoid checking back towards the central N atom
                if neighbor.GetAtomicNum() not in [1, 7, 8, 16, 15]: # Check for H and common heteroatoms
                    organic_neighbor = True
                elif neighbor.GetAtomicNum() == 6:
                     organic_neighbor = True # Also accept carbons
        if not organic_neighbor:
             return False, "Nitrogen is not bonded to three organic groups"


    return True, "Tertiary amine oxide detected"