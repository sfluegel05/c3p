"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone typically has two aromatic rings connected by a three-carbon α,β-unsaturated carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flexible SMARTS pattern for chalcones
    # Allows various substituent patterns on aromatic rings, connects to an α,β-unsaturated carbonyl
    aromatic_ring_smarts = "[$(c1ccccc1),$(c1ccccn1)]"  # Allowing heterocycles too
    enone_pattern = "[CX3](=O)[CX3]=[CX3]"  # General unsat linkage, α, β-unsat
    chalcone_smarts = f"{aromatic_ring_smarts}{enone_pattern}{aromatic_ring_smarts}"
    
    # Convert SMARTS string to pattern object
    chalcone_pattern = Chem.MolFromSmarts(chalcone_smarts)
    
    # Check for chalcone pattern in the molecule
    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Contains chalcone structure: two aromatic rings connected by α,β-unsaturated carbonyl"

    return False, "Does not contain chalcone structure: missing characteristic α,β-unsaturated carbonyl linkage"