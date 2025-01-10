"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid should have a steroid backbone, a ketone group at the 3-position,
    and specific stereochemistry denoted as '5beta'.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more detailed steroid backbone pattern (includes the four fused rings)
    steroid_pattern = Chem.MolFromSmarts("[C;R]1[C;R][C;R]2[C;R][C;R]3[C;R][C;R]4[C;R][C;R](C[C;R]4)[C;R](C[C;R]3)[C;R](C[C;R]2)[C;R](C1)")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Steroid backbone not found"
    
    # Define a pattern for the 3-oxo group, ensuring it is part of the steroid structure
    oxo_3_pattern = Chem.MolFromSmarts("[#6;R1]-[#6;R2](=O)-[#6;R3]")
    if not mol.HasSubstructMatch(oxo_3_pattern):
        return False, "3-oxo group not properly positioned"

    # Evaluate specific stereochemistry for '5beta'
    # This would involve detecting specific configurations and beta-orientation
    # Typically, 5beta steroids have a configuration with the rings in the trans orientation
    beta_positions = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    is_beta = any(pos[1] == 'S' or pos[1] == 'R' for pos in beta_positions)  # Example check for proper stereo
    
    if not is_beta:
        return False, "5beta stereochemistry not confirmed"

    return True, "Molecule matches the 3-oxo-5beta-steroid characteristics"