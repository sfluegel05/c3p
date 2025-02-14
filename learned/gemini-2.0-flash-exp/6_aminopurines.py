"""
Classifies: CHEBI:20706 6-aminopurines
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the core 6-aminopurine ring system
    core_ring_smarts = Chem.MolFromSmarts("c1nc2c(n1)ncn2")
    
     # Define the SMARTS pattern for the amino group at position 6
    amino_group_smarts = Chem.MolFromSmarts("[NX3;H0,H1,H2][#6]")

    if core_ring_smarts is None or amino_group_smarts is None:
        return False, "Invalid SMARTS string"

    # Combine the two SMARTS patterns
    combined_mol = Chem.CombineMols(core_ring_smarts, amino_group_smarts)

    # Check for the presence of the 6-aminopurine substructure
    if mol.HasSubstructMatch(combined_mol):
        return True, "Contains a 6-aminopurine moiety"
    else:
        return False, "Does not contain a 6-aminopurine moiety"