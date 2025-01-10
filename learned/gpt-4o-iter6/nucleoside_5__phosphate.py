"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose or deoxyribose sugar linked to a purine or pyrimidine base,
    with phosphorylation at the C-5 position of the sugar.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define generalized pattern for any ribose or deoxyribose containing nucleotide or deoxynucleotide
    sugar_pattern = Chem.MolFromSmarts("O[C@@]1([C@@H](O)C[C@@H](O)C1)[C@H]CO |C|")  # Uses generic C without specific stereocenters

    # Generalized phosphate group pattern that includes mono-, di-, tri-phosphate
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O |C|")  # Capture multiple conformations

    # General nucleotide base connection pattern
    base_patterns = [
        # Purines
        Chem.MolFromSmarts("n1cnc2c1[nH]cnc2"),  # Simplified purine base pattern
        # Pyrimidines like cytosine, uracil, thymine variations
        Chem.MolFromSmarts("n1ccn([C@@H]2O[C@H]3O[C@H](CO)[C@@H](O)[C@H]3O2)nc1=O")
    ]

    # Check for sugar, phosphate, and base pattern matches in the molecule
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "Does not contain a generalized ribose or deoxyribose sugar structure"

    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found attached to the sugar"

    if not any(mol.HasSubstructMatch(bp) for bp in base_patterns):
        return False, "No recognized nucleobase structure found linked via glycosidic bond"

    return True, "Contains a ribosyl or deoxyribosyl sugar, phosphate at 5', and a nucleobase"