"""
Classifies: CHEBI:26214 porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is characterized by a macrocyclic structure of four pyrrole nuclei 
    connected through the alpha-positions by four methine bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the porphyrin core
    porphyrin_pattern = Chem.MolFromSmarts('n1c(cc1)-c1nccc1-c1nccc1-c1nccc1')
    
    if porphyrin_pattern is None:
        return False, "Error in SMARTS pattern definition"

    # Check if the molecule matches the porphyrin pattern
    if mol.HasSubstructMatch(porphyrin_pattern):
        return True, "SMILES string contains a porphyrin core structure"
    else:
        return False, "No porphyrin core structure found in SMILES string"

# Example usage:
# print(is_porphyrins('C1=CC=C2C(=C1)C=CC=C2'))  # Example SMILES, not a porphyrin
# print(is_porphyrins('Oc1cccc(c1)-c1c2CCc(n2)c(-c2cccc(O)c2)c2ccc([nH]2)c(-c2cccc(O)c2)c2ccc1[nH]2'))  # Example SMILES for a porphyrin structure