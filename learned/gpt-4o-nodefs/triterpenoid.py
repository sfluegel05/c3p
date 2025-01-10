"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are typically composed of six isoprene units with a molecular formula
    often around C30H48 but can be highly functionalized.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight as it helps in balancing possible different substitutions
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 400 or mol_weight > 700:
        return False, f"Molecular weight of {mol_weight} is not typical for triterpenoids"

    # Check for the number of rings (generally 4 to 6)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4 or ring_count > 6:
        return False, f"Number of rings ({ring_count}) is not typical for triterpenoids"

    # Count carbon atoms to ensure appropriate isoprene skeleton
    num_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbon < 25 or num_carbon > 35:
        return False, f"Number of carbon atoms ({num_carbon}) not typical for triterpenoids"

    # Check for common triterpenoid skeletons (a simplified approximation using SMARTS)
    typical_triterpenoid_smarts = "C1CCC2C3C4C(C2)C(C1)C3C5C4C6C(C5)C(C6)=O"
    substruct_mol = Chem.MolFromSmarts(typical_triterpenoid_smarts)
    if not substruct_mol:
        return None, None  # Handle if the SMARTS pattern is incorrectly parsed
    
    backbone_match = mol.HasSubstructMatch(substruct_mol)
    
    if not backbone_match:
        return False, "Backbone structure doesn't match common triterpenoids"

    # Consider major functionalizations observed in triterpenoids
    functional_groups_smarts = ["[OH]", "[=O]", "[O;R]"]  # Hydroxyl, ketone, cyclic ethers
    for pattern in functional_groups_smarts:
        substruct_mol = Chem.MolFromSmarts(pattern)
        if not substruct_mol:
            return None, None  # Handle if the SMARTS pattern is incorrectly parsed

        if mol.HasSubstructMatch(substruct_mol):
            return True, "Structure matches characteristics of triterpenoids"
    
    return False, "No relevant functional groups detected"

# Example usage
# Test with a known triterpenoid SMILES
smiles = 'CC1(C)CCC2C3C(C(C4(C(C5(C(CC4)C(CCC=C(C)C)C)CC3)C)CC2)C1)C'
result = is_triterpenoid(smiles)
print(result)