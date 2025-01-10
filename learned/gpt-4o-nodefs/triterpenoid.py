"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 400 or mol_weight > 600:
        return False, f"Molecular weight of {mol_weight} is not typical for triterpenoids"

    # Check for the number of rings (generally 4 to 6)
    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4:
        return False, f"Number of rings ({ring_count}) is not enough for triterpenoids"

    # Count carbon atoms to ensure appropriate isoprene skeleton
    num_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbon < 25 or num_carbon > 35:
        return False, f"Number of carbon atoms ({num_carbon}) not typical for triterpenoids"

    # Check for common triterpenoid skeletons (a simplified approximation using SMARTS)
    # Example SMARTS for backbone structure -- needs customization for each case
    typical_triterpenoid_smarts = "[C;R]1[C;R][C;R][C;R]2[C;R][C;R]3[C;R]([C;R]4[C;R][C;R][C;R][C;R]12)[C;R]3"
    backbone_match = mol.HasSubstructMatch(Chem.MolFromSmarts(typical_triterpenoid_smarts))
    
    if not backbone_match:
        return False, "Backbone structure doesn't match common triterpenoids"

    # Consider major functionalizations observed in triterpenoids
    functional_groups_smarts = ["[OH]", "[=O]", "[O;R]", "[O;H0]"]  # Hydroxyl, ketone, cyclic ethers
    if any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in functional_groups_smarts):
        return True, "Structure matches characteristics of triterpenoids"
    else:
        return False, "No relevant functional groups detected"

# Example usage
# Test with a known triterpenoid SMILES
smiles = 'OC1C(C=2C(C3(C(C4(C(C(CC4)C(CCC=C(C)C)C)(CC3)C)C)CC2)C)CC1)(C)C'
result = is_triterpenoid(smiles)
print(result)