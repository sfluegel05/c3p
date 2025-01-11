"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Quinones have a fully conjugated cyclic dione structure, typically 
    derived from aromatic compounds, where an even number of -CH= groups 
    are replaced by -C(=O)- groups within an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a quinone structure
    quinone_pattern = Chem.MolFromSmarts("O=C1C=CC(=O)C=C1")
    
    # Check if the molecule matches the quinone pattern
    if mol.HasSubstructMatch(quinone_pattern):
        return True, "Molecule contains a fully conjugated cyclic dione structure indicative of a quinone"
    
    return False, "Molecule does not contain a quinone structural feature"