"""
Classifies: CHEBI:38958 indole alkaloid
"""
"""
Classifies: Indole Alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must contain an indole skeleton (a bicyclic structure consisting
    of a benzene ring fused to a pyrrole ring) and also contain functional groups
    or structural features characteristic of alkaloids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core indole skeleton SMARTS pattern including typical nitrogen configuration of indoles
    indole_pattern = Chem.MolFromSmarts('c1ccc2[nH]c([cH]c2c1)')
    
    # An additional SMARTS pattern that includes a broad definition of nitrogen functionalities in alkaloids
    # This assumes indole structures with additional nitrogen like substituents common in alkaloids
    indole_alkaloid_pattern = Chem.MolFromSmarts('c1ccc2[nR]c([aH])[nR][aH]c2c1')

    # Match using both primary pattern and the broadened indole alkaloid pattern
    if mol.HasSubstructMatch(indole_pattern) or mol.HasSubstructMatch(indole_alkaloid_pattern):
        
        # Further check for the presence of core alkaloid features
        nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
        if len(nitrogen_atoms) > 1:  # Expect more than one nitrogen in more complex structures
            return True, "Contains indole skeleton with additional alkaloid features"

        return False, "Contains indole but lacks additional features typical of alkaloids"
    else:
        return False, "No indole skeleton found"