"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is characterized by the presence of multiple amine groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for multiple amine groups (-NH2, -NH-)
    amine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")  # Simple primary/secondary amine pattern
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    
    # Check if there are at least two distinct amine groups
    if len(amine_matches) < 2:
        return False, f"Found {len(amine_matches)} amine group(s), need at least 2"

    # Evaluate the distribution of amine groups across the structure
    # In a polyamine, amine groups can be spaced across aliphatic or aromatic chains
    atom_indices = [match[0] for match in amine_matches]
    
    # Check if they are potentially connected by analyzing distances
    # Placeholder logic to evaluate or visualize sequence - actual pattern might include assessing chain separation
    # For simplicity, assume that multiple distinct amines are enough without further spatial analysis

    return True, f"Contains multiple amine groups: {len(amine_matches)} found"

# The function is_polyamine can now be used to classify polyamines based on the SMILES input.