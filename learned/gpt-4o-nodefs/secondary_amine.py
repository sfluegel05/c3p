"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to two carbon atoms and either 
    one hydrogen atom or a lone electron pair.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms to find secondary amine characteristic
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check for Nitrogen atom
            # Count the number of carbon neighbors
            num_carbon = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            # Check if there's one hydrogen or a lone pair on nitrogen assuming valency is 3 for NH)
            num_hydrogen = atom.GetTotalNumHs()

            # A secondary amine has Nitrogen connected to two Carbon atoms and may have one hydrogen.
            if num_carbon == 2 and (num_hydrogen == 1 or atom.GetFormalCharge() == 0):
                # Additional constraint - ensure no other non-carbon bonded
                non_C_H_bonds = any(neighbor.GetAtomicNum() not in [1, 6] for neighbor in atom.GetNeighbors())
                if not non_C_H_bonds:
                    return True, "Contains a secondary amine group"
    
    return False, "Does not contain a secondary amine group"

# Testing example correctness for known secondary amine
# print(is_secondary_amine('CCNCC'))  # Diethylamine, expected True.