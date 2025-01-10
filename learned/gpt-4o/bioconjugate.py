"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate consists of at least two biological molecules covalently linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Simple fragments for known bioconjugates
    patterns = {
        "peptide_bond": Chem.MolFromSmarts("N[C@H](C)C(=O)"),
        "disulfide_bond": Chem.MolFromSmarts("S-S"),
        "glutathione_like": Chem.MolFromSmarts("N[C@H](CC(=O)NCCS)C(=O)NCC(=O)O"),
        "coenzyme_a_link": Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)C@[O,P]")
    }
    
    # Flag to check if we match any bioconjugate-like pattern
    matched_patterns = []
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            matched_patterns.append(name)

    # We check for the presence of at least two distinct biological units
    if len(matched_patterns) >= 2:
        return True, f"Contains patterns: {', '.join(matched_patterns)}"

    # For failed cases, provide the reason for classification based on detected patterns
    if matched_patterns:
        return False, f"Partially matched patterns, only found: {', '.join(matched_patterns)}"

    # More complex logic might be needed to cover diverse examples encountered in real applications
    return False, "No definitive bioconjugate patterns found"