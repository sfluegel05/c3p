"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids contain an aromatic ring connected to a propane chain 
    or similar structures like flavonoids, coumarins, etc.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define aromatic systems that can be part of phenylpropanoids.
    aromatic_system_patterns = [
        Chem.MolFromSmarts("c1ccccc1"),  # Simple benzene
        Chem.MolFromSmarts("c1ccc2ccccc2c1"),  # Naphthalene or similar
        Chem.MolFromSmarts("c1cc(c2ccccc2)ccc1"),  # Phenyl derivative
    ]
    
    # Check for presence of these aromatic systems
    has_aromatic_system = any(mol.HasSubstructMatch(pattern) for pattern in aromatic_system_patterns)
    
    if not has_aromatic_system:
        return False, "No aromatic system characteristic of phenylpropanoids found"
    
    # Look for potential phenylpropanoid backbone patterns
    phenylpropanoid_backbone_patterns = [
        Chem.MolFromSmarts("c-C-C-C"),  # Simple phenylpropane backbone
        Chem.MolFromSmarts("c-CO"),  # Coumarin-related 
        Chem.MolFromSmarts("c-O-c"),  # Flavonoid
    ]
    
    # Check for presence of these patterns
    has_phenylpropanoid_backbone = any(mol.HasSubstructMatch(pattern) for pattern in phenylpropanoid_backbone_patterns)
    
    if not has_phenylpropanoid_backbone:
        return False, "No phenylpropanoid backbone found"
    
    # Check for functional groups common in phenylpropanoids
    functional_group_patterns = [
        Chem.MolFromSmarts("c-[OH]"),  # Phenolic OH group
        Chem.MolFromSmarts("c-CO"),    # Carbonyl
        Chem.MolFromSmarts("c-COC"),   # Methoxy
        Chem.MolFromSmarts("c=O"),     # Part of the carbonyl system in coumarins
    ]
    
    has_functional_group = any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns)
    
    if not has_functional_group:
        return False, "No common phenylpropanoid functional groups found"

    return True, "Structure contains aromatic system characteristic of a phenylpropanoid and necessary functional groups."