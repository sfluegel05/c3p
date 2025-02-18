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

    # Define the core penam structure with key substituents
    penam_core = Chem.MolFromSmarts("[*]C1([CH3])[CH3][S][C@@H]2[C@@H](NC(=O)[*])[C@](=O)(O)N12")
    if not mol.HasSubstructMatch(penam_core):
        return False, "Missing penam core with two methyl groups and carboxylate"

    # Check for carboxamido group (R-N-C(=O)-R') attached to the beta-lactam nitrogen
    carboxamido = Chem.MolFromSmarts("[NX3][C](=O)")
    matches = mol.GetSubstructMatches(carboxamido)
    if not matches:
        return False, "No carboxamido group found"

    # Ensure the carboxamido is attached to the beta-lactam nitrogen in the core
    core_matches = mol.GetSubstructMatches(penam_core)
    for cm in core_matches:
        beta_lactam_n_idx = cm[4]  # Index of the beta-lactam nitrogen in the core match
        for (n_idx, c_idx) in matches:
            if n_idx == beta_lactam_n_idx:
                return True, "Penicillin core with carboxamido group found"

    return False, "Carboxamido group not attached to beta-lactam nitrogen"