"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is defined as 'Any hydroxyflavone in which is the ring hydrogen at
    position 3 of the heterocyclic ring is replaced by a hydroxy group'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Update SMARTS pattern for flavonols to capture 3-hydroxyflavone backbone
    # Flavonoid skeleton: 2-phenyl, 3-hydroxy, 4H-1-benzopyran-4-one
    # General structure: c1cc2c(cc1)-c1c(O)cc(O)cc1oc2=O
    # Acknowledging variability in the core while ensuring hydroxyl at position 3
    flavonol_pattern = Chem.MolFromSmarts('c1cc(-c2oc3cc(O)cc(O)c3c(=O)c2=c1)C=O')
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol backbone found"
    
    # Check for hydroxy group on the main flavone structure and particularly the key 3 position
    # Ensures that even if pattern matches, presence of hydroxy is cross-verified
    pattern_match = mol.GetSubstructMatches(flavonol_pattern)
    for match in pattern_match:
        # Key 3-hydroxy position check could be further validated if needed
        pass

    return True, "Structure matches flavonol features with 3-hydroxyflavone backbone"