"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound contains at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop through atoms to check for fluorine with a carbon neighbor
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 9:  # Fluorine atomic number is 9
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):  # Carbon atomic number is 6
                return True, "Molecule contains at least one carbon-fluorine bond"
    
    return False, "No carbon-fluorine bond found in the molecule"

# Example usage:
# result, reason = is_organofluorine_compound("C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=C(C=C3)F)O[C@H]1CN(C)C(=O)C)[C@@H](C)CO")
# print(result, reason)