"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds, often forming complex glucosides of catechol and pyrogallol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for complex polyphenolic structures, including catechol and pyrogallol motifs
    polyphenol_pattern = Chem.MolFromSmarts("c1(c(O)c(O)c(O)c1)~c2cc(O)c(O)c(O)c2")  # Larger patterns for complex phenolic structures
    polyphenol_matches = mol.GetSubstructMatches(polyphenol_pattern)
    if len(polyphenol_matches) < 1:
        return False, "Not enough complex polyphenolic groups detected"

    # Verify presence of sufficient hydroxyl groups indicative of polyphenolic nature
    o_count = sum(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() > 0 for atom in mol.GetAtoms())
    if o_count < 8:  # Increase threshold to capture extensive hydroxylation
        return False, f"Too few hydroxyl groups, found {o_count}"

    # Optional check for glucosidic linkages (e.g., -O-C-C-O-)
    glucoside_pattern = Chem.MolFromSmarts("O[C@H|@@H]([C@H|@@H]1O)[C@H|@@H](O1)C")  # A more generalized glucoside pattern
    if not mol.HasSubstructMatch(glucoside_pattern):
        return False, "No glucosidic linkage detected"

    return True, "Polyphenolic structure with catechol/pyrogallol motifs and glucosidic linkages detected"