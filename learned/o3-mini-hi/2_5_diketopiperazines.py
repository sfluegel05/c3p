"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines (Any piperazinone that has a piperazine-2,5-dione skeleton.)
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is defined as any piperazinone that has a piperazine-2,5-dione
    skeleton. The core structure contains a six-membered ring composed of two nitrogen atoms
    and two carbonyl groups located at positions 2 and 5 of the ring.
    
    The simplest representation of this scaffold is: O=C1NC(=O)CN1

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the piperazine-2,5-dione core.
    # This pattern looks for a six-membered ring with two amide (=O) groups attached.
    diketopiperazine_smarts = "O=C1NC(=O)CN1"
    pattern = Chem.MolFromSmarts(diketopiperazine_smarts)
    if pattern is None:
        return False, "Error in SMARTS pattern definition"
        
    # Check if the molecule contains the diketopiperazine substructure
    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains a piperazine-2,5-dione (2,5-diketopiperazine) skeleton"
    else:
        return False, "Molecule does not contain the required piperazine-2,5-dione skeleton"
        
# Example usage (uncomment to test):
# smiles_example = "O=C1NC(=O)CN1"  # piperazine-2,5-dione, simplest case.
# result, reason = is_2_5_diketopiperazines(smiles_example)
# print(result, reason)