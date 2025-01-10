"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are complex polyphenolic compounds forming glucosides of catechol and pyrogallol.

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

    # Look for catechol and pyrogallol or similar aromatic polyphenolic subunits
    polyphenol_patterns = [
        Chem.MolFromSmarts("c1(c(O)c(O)c(O)c1)~"),   # General pattern to match polyphenol ethers
        Chem.MolFromSmarts("c1cc(O)c(O)c(c1)O"),  # Catechol
        Chem.MolFromSmarts("c1c(O)c(O)cc(O)c1"),  # Pyrogallol
        Chem.MolFromSmarts("c1cc(c(c(c1O)O)O)O")  # Extended pyrogallol
    ]

    # Check each pattern to find if any complex polyphenolic structure matches
    found_polyphenol = any(mol.HasSubstructMatch(pattern) for pattern in polyphenol_patterns if pattern)

    if not found_polyphenol:
        return False, "No complex polyphenolic groups detected; exhibits simpler phenolic structure"

    # Count number of hydroxyl and ether groups
    o_count = sum(atom.GetAtomicNum() == 8 for atom in mol.GetAtoms())
    if o_count < 10:  # Tannins are rich in oxygen due to their polyhydroxy nature
        return False, f"Too few hydroxyl/ether groups detected, found {o_count}"

    # Look for glucoside linkages or ester bonds typical of tannins
    glucoside_pattern = Chem.MolFromSmarts("O[C@H1]([C@H1](O)[C@H1](O)[C@H1](O)")[C@H1](O)[C@H1](O)C)~") # Example glucose substructure
    found_glucoside = mol.HasSubstructMatch(glucoside_pattern)

    if not found_glucoside:
        return False, "No glucosidic linkages detected"

    return True, "Complex polyphenolic structure with glucosidic linkages and catechol/pyrogallol motifs detected"