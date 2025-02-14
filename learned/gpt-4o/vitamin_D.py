"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are identified by their hydroxy seco-steroid backbone
    and specific structural features that are characteristic of this family of molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved seco-steroid pattern for vitamin D
    # This pattern looks for the characteristic cleavage of the B-ring in the steroid nucleus
    seco_steroid_pattern = Chem.MolFromSmarts("C1C=C2C[C@H](C[C@@H](O)C2)O1") # Improved pattern for broken B-ring
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No characteristic seco-steroid structure found"

    # Ensure the presence of necessary hydroxyl groups in typical positions for vitamin D
    hydroxyl_groups_required = [
        "[C@H](O)",  # A chiral center with hydroxyl, indicative of position-specific hydroxylation
    ]
    hydroxyl_match_count = 0
    for pattern in hydroxyl_groups_required:
        group_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(group_pattern):
            hydroxyl_match_count += 1
    
    if hydroxyl_match_count < 2:
        return False, f"Found {hydroxyl_match_count} hydroxyl groups, characteristic requires more"

    # Check for the existence of long alkyl side chain and stereochemistry typical of vitamin D
    alkyl_chain_pattern = Chem.MolFromSmarts("C(C)CCCC(C)")
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "Missing typical long alkyl side chain"

    # If all characteristics are matched, classify as vitamin D
    return True, "Molecule matches key structural features of vitamin D"