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
    
    # Define specific patterns to recognize components within bioconjugates
    patterns = {
        "peptide_bond": Chem.MolFromSmarts("N[C@H](C)C(=O)"),  # Typical amide in peptides
        "thioester_bond": Chem.MolFromSmarts("C(=O)S"),  # Common in acyl-CoA
        "disulfide_bond": Chem.MolFromSmarts("S-S"),  # Disulfide bonds in proteins
        "glutathione_motif": Chem.MolFromSmarts("N[C@H](CC(=O)NCCS)C(=O)NCC(=O)")  # Part of glutathione
    }

    # Other potential patterns that may hint at biological functionality
    accessory_patterns = {
        "phosphate_group": Chem.MolFromSmarts("P(=O)(O)O"),  # Nucleotide or ATP-like
        "heterocyclic_nitrogen": Chem.MolFromSmarts("n"),  # Purine/pyrimidine-like
        "ether_linkage": Chem.MolFromSmarts("C-O-C"),  # Often seen in glycolipids
        "sulfur_containing": Chem.MolFromSmarts("S")  # Can indicate methionine, cysteine, CoA
    }
    
    # Track which primary patterns we match
    matched_patterns = set()
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            matched_patterns.add(name)
    
    # Track which accessory patterns we match
    accessory_matches = set()
    for name, pattern in accessory_patterns.items():
        if mol.HasSubstructMatch(pattern):
            accessory_matches.add(name)

    # Determine if there are at least two distinctive substructures suggesting bioconjugation
    if len(matched_patterns) >= 2 or (len(matched_patterns) == 1 and len(accessory_matches) > 1):
        return True, f"Contains patterns: {', '.join(matched_patterns.union(accessory_matches))}"

    # Provide reasoning for non-bioconjugate classification
    if matched_patterns or accessory_matches:
        return False, f"Partially matched patterns, found: {', '.join(matched_patterns.union(accessory_matches))}"

    # Return explanation for no pattern matches
    return False, "No definitive bioconjugate patterns found"