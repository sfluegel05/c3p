"""
Classifies: CHEBI:23437 cyanohydrin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin has a hydroxyl and a nitrile group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the cyanohydrin pattern: a carbon bonded to -OH and -C#N
    cyanohydrin_pattern = Chem.MolFromSmarts("[C;D3]([OH])([#7]#[#6])")
    if mol.HasSubstructMatch(cyanohydrin_pattern):
        return True, "Contains cyanohydrin moiety (-C(OH)(C#N)-)"
    else:
        return False, "Does not contain cyanohydrin moiety"

# Example usage:
# print(is_cyanohydrin("O[C@@H](C#N)c1ccccc1"))  # Should return (True, "...Contains cyanohydrin moiety...")