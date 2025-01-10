"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are typically characterized by long conjugated polyene chains,
    presence of cyclic or acyclic end groups, and various functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS patterns for key features of carotenoids
    long_conjugated_chain_pattern = Chem.MolFromSmarts("C(=C-C=C){4,}")  # Conjugated chains
    cyclic_end_pattern = Chem.MolFromSmarts("C1=CC=CC=C1")  # Simple aromatic cycle example
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    epoxide_pattern = Chem.MolFromSmarts("[C]-[O]-[C]")
    
    # Check for polyene chain
    if mol.HasSubstructMatch(long_conjugated_chain_pattern):
        # Further check for specific functional groups or cyclic structures
        has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
        has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)
        has_epoxide = mol.HasSubstructMatch(epoxide_pattern)
        has_cycle = mol.HasSubstructMatch(cyclic_end_pattern)
        
        # Compile characteristics
        functional_groups = has_hydroxyl or has_carbonyl or has_epoxide
        if functional_groups or has_cycle:
            return True, "Contains long conjugated chain and functional/cyclic end groups typical of carotenoids"
    
    return False, "Does not match carotenoid structural features or complexity"

# Example usage with one of the given carotenoid structures
example_smiles = "CC(\C=C\C=C(C)\C=C\C=C(C)\C=C\[C@@H](O)C(C)(C)C)=C/C=C/C=C(C)\C=C/C=C(C)/C=C/C1=C(C)C(O)[C@@H](O)CC1(C)C"
print(is_carotenoid(example_smiles))