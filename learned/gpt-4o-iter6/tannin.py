"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds often forming complex glucosides of catechol and pyrogallol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polyphenolic structure: multiple aromatic rings with hydroxyls
    # Relax the pattern to identify various polyphenolic forms
    polyphenol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1")  # Allow rings with multiple hydroxyl groups
    polyphenol_matches = mol.GetSubstructMatches(polyphenol_pattern)
    if len(polyphenol_matches) < 1:
        return False, "Not enough polyphenolic groups detected"

    # Count hydroxyl groups to establish polyphenolic nature
    o_count = sum(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0 for atom in mol.GetAtoms())
    if o_count < 6:
        return False, f"Too few hydroxyl groups, found {o_count}"

    return True, "Polyphenolic structure with sufficient hydroxyl groups detected"