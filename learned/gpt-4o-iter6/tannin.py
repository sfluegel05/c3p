"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are complex polyphenolic compounds often forming glucosides of catechol and pyrogallol.

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

    # Expanded polyphenolic substructure matching
    # Detecting various catechol and pyrogallol-like motifs
    polyphenol_patterns = [
        Chem.MolFromSmarts("c1(c(O)c(O)c(O)c1)~c2cc(O)c(O)c(O)c2"),  # Extended catechol
        Chem.MolFromSmarts("c1cc(O)c(O)c(c1)O"),  # Catechol
        Chem.MolFromSmarts("c1c(O)c(O)cc(O)c1")  # Pyrogallol
    ]
    
    # Check for presence of complex polyphenolic groups
    found_polyphenol = False
    for pattern in polyphenol_patterns:
        if mol.HasSubstructMatch(pattern):
            found_polyphenol = True
            break
    if not found_polyphenol:
        return False, "No complex polyphenolic groups detected; exhibits simpler phenolic structure"

    # Verify presence of sufficient hydroxyl groups
    o_count = sum(atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 for atom in mol.GetAtoms())
    if o_count < 15:  # Modified threshold to capture extensive hydroxylation typical of tannins
        return False, f"Too few hydroxyl groups detected, found {o_count}"

    # Expanded search for glucosidic linkages and similar motifs
    glucoside_patterns = [
        Chem.MolFromSmarts("O[C@H|@@H]([C@H|@@H]1O)[C@H|@@H](O1)"),  # Primary glucoside pattern
        Chem.MolFromSmarts("O=C(O[C@H|@@H]1OC[C@H|@@H](O)[C@H|@@H](O)[C@H|@@H](O)[C@H|@@H]1)") # Esterified glucose
    ]
    found_glucoside = False
    for gluc_pattern in glucoside_patterns:
        if mol.HasSubstructMatch(gluc_pattern):
            found_glucoside = True
            break
    if not found_glucoside:
        return False, "No glucosidic linkage or related ester motifs detected"

    return True, "Complex polyphenolic structure with glucosidic linkages and catechol/pyrogallol motifs detected"