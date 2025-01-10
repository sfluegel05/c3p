"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is defined as 'Any hydroxyflavone in which the ring hydrogen at
    position 3 of the heterocyclic ring is replaced by a hydroxy group'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavonol, False otherwise
        str: Reason for the classification
    """
    
    # Parse SMILES to a RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the Smart pattern for a flavonol structure.
    # It includes specific position for 3-hydroxy group in the flavone backbone.
    flavonol_pattern = Chem.MolFromSmarts('O=c1oc2cc(O)cc(O)c2c(O)c1')

    # Check if the structure contains a flavonol pattern
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol backbone detected with 3-hydroxy group"

    return True, "Structure matches 3-hydroxyflavone backbone consistent with flavonol"