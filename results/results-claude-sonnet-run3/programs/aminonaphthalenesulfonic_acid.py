from rdkit import Chem
from rdkit.Chem import AllChem

def is_aminonaphthalenesulfonic_acid(smiles: str):
    """
    Determines if a molecule is an aminonaphthalenesulfonic acid, defined as
    a naphthalenesulfonic acid having at least one amino substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aminonaphthalenesulfonic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for naphthalene core
    naphthalene_pattern = Chem.MolFromSmarts('c1cccc2ccccc12')
    if not mol.HasSubstructMatch(naphthalene_pattern):
        return False, "No naphthalene core found"

    # Check for sulfonic acid group
    sulfonic_acid_pattern = Chem.MolFromSmarts('S(=O)(=O)O')
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group found"

    # Check for amino group (both primary and secondary amines)
    primary_amine_pattern = Chem.MolFromSmarts('[NH2]')
    secondary_amine_pattern = Chem.MolFromSmarts('[NH;!H0]')
    
    has_primary_amine = mol.HasSubstructMatch(primary_amine_pattern)
    has_secondary_amine = mol.HasSubstructMatch(secondary_amine_pattern)

    if not (has_primary_amine or has_secondary_amine):
        return False, "No amino group found"

    # Count substituents
    num_sulfonic = len(mol.GetSubstructMatches(sulfonic_acid_pattern))
    num_primary_amino = len(mol.GetSubstructMatches(primary_amine_pattern))
    num_secondary_amino = len(mol.GetSubstructMatches(secondary_amine_pattern))

    reason = f"Found naphthalene core with {num_sulfonic} sulfonic acid group(s) and "
    if num_primary_amino > 0 and num_secondary_amino > 0:
        reason += f"{num_primary_amino} primary and {num_secondary_amino} secondary amino group(s)"
    elif num_primary_amino > 0:
        reason += f"{num_primary_amino} primary amino group(s)"
    else:
        reason += f"{num_secondary_amino} secondary amino group(s)"

    return True, reason
# Pr=None
# Recall=None