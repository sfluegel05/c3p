"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: Penicillin (CHEBI:17334)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    Penicillin has a penam core with two methyl groups at position 2,
    a carboxylate at position 3, and a carboxamido group at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for two methyl groups attached to sulfur
    methyl_sulfur = Chem.MolFromSmarts("[S]C(C)(C)")
    if not mol.HasSubstructMatch(methyl_sulfur):
        return False, "No sulfur with two methyl groups"

    # Check adjacent carbon has carboxylate (C(=O)O)
    carboxylate = Chem.MolFromSmarts("[S]C(C)(C)[C](=O)[OX1]")
    if not mol.HasSubstructMatch(carboxylate):
        return False, "No carboxylate adjacent to sulfur"

    # Check for beta-lactam ring (4-membered ring with N-C=O)
    beta_lactam = Chem.MolFromSmarts("[n]1[c](=O)[c][c]1")
    if not mol.HasSubstructMatch(beta_lactam):
        beta_lactam_alt = Chem.MolFromSmarts("[NX3]1C(=O)CC1")  # Alternative pattern
        if not mol.HasSubstructMatch(beta_lactam_alt):
            return False, "No beta-lactam ring"

    # Find the beta-lactam nitrogen and check for attached carboxamido group
    carboxamido = Chem.MolFromSmarts("[NX3][C](=O)")
    beta_lactam_n = Chem.MolFromSmarts("[NX3]C(=O)")
    n_matches = mol.GetSubstructMatches(beta_lactam_n)
    if not n_matches:
        return False, "No beta-lactam nitrogen found"

    # Check if any beta-lactam nitrogen has a carboxamido substituent
    for n_idx, c_idx in mol.GetSubstructMatches(carboxamido):
        for n_match in n_matches:
            if n_idx == n_match[0]:
                return True, "Penicillin core with required substituents found"

    return False, "Carboxamido group not attached to beta-lactam nitrogen"