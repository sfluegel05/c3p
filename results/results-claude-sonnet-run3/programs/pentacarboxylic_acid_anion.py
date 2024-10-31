from rdkit import Chem
from rdkit.Chem import AllChem

def is_pentacarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a pentacarboxylic acid anion.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pentacarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Count number of carboxylate anions (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts('[C](=[O])[O-]')
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    num_carboxylate = len(carboxylate_matches)

    # Count number of carboxylic acids (-COOH)
    carboxyl_pattern = Chem.MolFromSmarts('[C](=[O])[OH]')
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    num_carboxyl = len(carboxyl_matches)

    total_carboxy = num_carboxylate + num_carboxyl

    # Must have exactly 5 total carboxy groups (sum of -COO- and -COOH)
    if total_carboxy != 5:
        return False, f"Found {total_carboxy} carboxy groups, need exactly 5"

    # Must have at least 1 carboxylate anion
    if num_carboxylate == 0:
        return False, "No carboxylate anions found"

    return True, f"Contains {num_carboxylate} carboxylate anions and {num_carboxyl} carboxylic acid groups"
# Pr=None
# Recall=None