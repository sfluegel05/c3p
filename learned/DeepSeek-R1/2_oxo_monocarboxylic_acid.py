"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:36586 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid has exactly one carboxylic acid group with a ketone on the adjacent carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(carboxylic_matches) != 1:
        return False, f"Found {len(carboxylic_matches)} carboxylic acid groups, need exactly 1"

    # Get the carbon atom of the carboxylic acid group
    carb_carbon_idx = carboxylic_matches[0][0]
    carb_carbon = mol.GetAtomWithIdx(carb_carbon_idx)

    # Check adjacent carbons for a ketone group (C=O)
    found_oxo = False
    for neighbor in carb_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Check only carbon neighbors
            # Look for double bond to oxygen (ketone)
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetAtomicNum() == 8:  # Oxygen atom in the double bond
                        found_oxo = True
                        break
            if found_oxo:
                break

    if not found_oxo:
        return False, "No oxo group adjacent to carboxylic acid carbon"

    return True, "Monocarboxylic acid with 2-oxo substituent"