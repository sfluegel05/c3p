"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Identifies conjugated polyene chains and typical carotenoid functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for long conjugated chain
    conjugated_chain_pattern = Chem.MolFromSmarts("C=C(-*):*-*...")
    if not mol.HasSubstructMatch(conjugated_chain_pattern):
        return False, "No long conjugated chain detected"

    # Check length requirement - minimal number of conjugated double bonds
    num_conjugated_dbs = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2.0)
    if num_conjugated_dbs < 9:
        return False, f"Conjugated system too short; found {num_conjugated_dbs} double bonds"

    # Check for typical carotenoid functional groups (like hydroxyl, carbonyl)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    carbonyl_pattern = Chem.MolFromSmarts("[C]=[O]")
    
    if mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(carbonyl_pattern):
        return True, "Contains long conjugated system and typical carotenoid functional groups"

    # Attempt to identify rings or end groups if possible
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return True, "Contains long conjugated system with ring structures"

    return True, "Contains long conjugated system characteristic of carotenoids"


# Example usage
# Note: Replace the SMILES string below with actual example SMILES for testing
example_smiles = "CC(\C=C\C=C(C)\C=C\c1c(C)ccc(C)c1C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/c1c(C)ccc(C)c1C"
print(is_carotenoid(example_smiles))