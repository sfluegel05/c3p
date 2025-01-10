"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are typically composed of six isoprene units with a molecular formula
    close to C30H48, often possessing multiple rings.

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

    # Check molecular formula
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    if not formula.startswith('C30H'):
        return False, f"Molecular formula is {formula}, not consistent with common triterpenoids"

    # Check for the number of rings
    ring_count = Chem.rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4 or ring_count > 6:
        return False, f"Number of rings ({ring_count}) is not typical for triterpenoids"

    # It is difficult to determine isoprene units directly from the SMILES,
    # but checking for the number of carbon atoms can be a proxy check
    num_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (25 <= num_carbon <= 35):
        return False, f"Number of carbon atoms ({num_carbon}) not typical for triterpenoids"

    # Check for characteristic functional groups
    potential_functional_groups = ['O', 'OH', '=O', 'OC', 'CO']  # Sparse check for simple functional groups
    relevant_functional_groups_present = any(mol.HasSubstructMatch(Chem.MolFromSmarts(group)) for group in potential_functional_groups)

    if not relevant_functional_groups_present:
        return False, "No relevant functional groups detected"

    return True, "Structure matches characteristics of triterpenoids"

# Example usage
# Test with a known triterpenoid SMILES
smiles = 'OC1C(C=2C(C3(C(C4(C(C(CC4)C(CCC=C(C)C)C)(CC3)C)C)CC2)C)CC1)(C)C'
result = is_triterpenoid(smiles)
print(result)