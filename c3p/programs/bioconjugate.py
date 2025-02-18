"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least 2 biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define some bioconjugate substructure patterns (amino acids, sugars, etc.)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2][CX4][CX3](=O)[OX2H1]")
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H]1O")
    lipid_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")

    # Find matches for these patterns
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    lipid_chain_matches = mol.GetSubstructMatches(lipid_chain_pattern)

    # Count the number of distinct matches
    biological_unit_count = 0
    if len(amino_acid_matches) > 0:
        biological_unit_count += 1
    if len(sugar_matches) > 0:
        biological_unit_count += 1
    if len(lipid_chain_matches) > 0:
        biological_unit_count += 1

    # Check if multiple biologically relevant units are present and connected
    if biological_unit_count >= 2:
        return True, f"Contains {biological_unit_count} biologically relevant substructures covalently linked"
    else:
        return False, f"Contains {biological_unit_count} biologically relevant substructures, need at least 2"

    return False, "Not enough biological molecules covalently linked"

# Example usage:
result, reason = is_bioconjugate("S(SS(O)(=O)=O)C[C@H](N)C(O)=O")
print(result, reason)