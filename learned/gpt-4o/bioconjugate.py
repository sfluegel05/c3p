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
    
    # Expanded patterns for detecting common bioconjugate motifs
    patterns = {
        "peptide_bond": Chem.MolFromSmarts("N[C@H](C)C(=O)"),
        "disulfide_bond": Chem.MolFromSmarts("S-S"),
        "glutathione_like": Chem.MolFromSmarts("N[C@H](CC(=O)NCCS)C(=O)NCC(=O)O"),
        "coenzyme_a_linkage": Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)"),
        "thioester_bond": Chem.MolFromSmarts("C(=O)S"),
        "ester_bond": Chem.MolFromSmarts("C(=O)O"),
        "amide_bond": Chem.MolFromSmarts("C(=O)N"),
        "heterocyclic_nitrogen": Chem.MolFromSmarts("n")
    }
    
    # Track which patterns we match
    matched_patterns = set()
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            matched_patterns.add(name)

    # Determine if there are at least two distinctive substructures suggesting bioconjugation
    if len(matched_patterns) >= 2:
        return True, f"Contains patterns: {', '.join(matched_patterns)}"

    # For failed cases, provide reasoning based on detected patterns
    if matched_patterns:
        return False, f"Partially matched patterns, only found: {', '.join(matched_patterns)}"

    # Provide explanation when no patterns are matched
    return False, "No definitive bioconjugate patterns found"