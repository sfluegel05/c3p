"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    Steroids have a cyclopenta[a]phenanthrene skeleton (three fused six-membered rings and one five-membered ring),
    with methyl groups at C-10 and C-13, and an optional alkyl group at C-17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid core structure using SMARTS pattern
    # This pattern approximates the fused four-ring system (cyclopenta[a]phenanthrene)
    steroid_core = Chem.MolFromSmarts("[C@]12[C@@H]([C@@H]3[C@H](CC1)C[C@H]([C@]4([C@H]([C@@H]23)CCCC4)C)C)C")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Missing cyclopenta[a]phenanthrene core"

    # Check for at least two methyl groups (C-10 and C-13)
    methyl_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CH3]"))
    if len(methyl_matches) < 2:
        return False, f"Found only {len(methyl_matches)} methyl groups, need at least 2"

    # Check for possible alkyl group at C-17 (any non-methyl carbon chain)
    c17_alkyl = mol.GetSubstructMatches(Chem.MolFromSmarts("[C][C]"))  # Simple check for at least a two-carbon chain
    if not c17_alkyl:
        return False, "No alkyl group detected at C-17"

    return True, "Contains steroid core with methyl groups and C-17 alkyl substituent"