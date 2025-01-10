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
    
    # Porphyrin rings have four nitrogen atoms in a cyclic, macrocyclic structure
    # The updated SMARTS accounts for the macrocyclic porphyrin with pyrrole units and connections
    # Patterns: c1ccc2nccc2c1 - 4 times in a macrocycle
    porphyrin_motif_pattern = Chem.MolFromSmarts('[n]1ccc2[n]ccc2c3[n]ccc4[n]ccc1c34')
    if porphyrin_motif_pattern is None:
        return False, "Error in SMARTS pattern definition"

    # Check if the molecule matches the porphyrin pattern
    if mol.HasSubstructMatch(porphyrin_motif_pattern):
        return True, "SMILES string contains a porphyrin core structure with cyclic pyrrole connections"
    else:
        return False, "No porphyrin core structure found in SMILES string"

# Example usage:
# print(is_porphyrins('Oc1cccc(c1)-c1c2CCc(n2)c(-c2cccc(O)c2)c2ccc([nH]2)c(-c2cccc(O)c2)c2ccc1[nH]2'))  # Example SMILES for a porphyrin structure