"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids typically have long conjugated polyene chains and may include 
    functional groups such as hydroxyl or carbonyl groups.

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
    
    # Define SMARTS patterns for characteristics of carotenoids
    conjugated_chain_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=[C]-[C]=[C]-[C]=[C]")  # at least three conjugated double bonds as minimum example
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")

    # Check for conjugated chain
    if not mol.HasSubstructMatch(conjugated_chain_pattern):
        return False, "No long conjugated polyene chain detected; fundamental for carotenoids"

    # Check for functional groups
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_carbonyl = mol.HasSubstructMatch(carbonyl_pattern)

    if has_hydroxyl or has_carbonyl:
        return True, "Contains long conjugated polyene chain and carotenoid-like functional groups"
    
    # Check for generic carotenoid structure based on size and molecular weight (as carotenoids are typically large)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)  # Large number implies complexity typical of carotenoids
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Arbitrary cutoff values based on typical carotenoids' properties
    if n_rotatable > 20 and mol_wt > 500:
        return True, "Contains long conjugated polyene chain typical of carotenoids"

    return False, "Lacks typical carotenoid functional groups and size"

# Example usage (replace with actual example SMILES for testing)
example_smiles = "CC(=C)CCCC(=C)C(=C)C(C)=CC"
print(is_carotenoid(example_smiles))