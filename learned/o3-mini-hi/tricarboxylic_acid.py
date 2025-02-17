"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: Tricarboxylic Acid (An oxoacid containing three carboxy groups.)
 
This improved version uses two SMARTS patterns to detect carboxyl groups 
in either protonated ([CX3](=O)[OX2H]) or ionized ([CX3](=O)[O-]) form and 
then checks that exactly three unique carboxyl carbons are detected.
It also excludes molecules that contain metals or other non‐organic atoms,
since many false positives were found from salts or complexes.
Note that this classification is heuristic.
"""

from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is defined as an oxoacid containing three carboxy groups.
    This function checks if the input molecule has exactly three carboxyl groups
    (in protonated or deprotonated form) and that it is composed solely of
    typical non-metal organic atoms (C, H, N, O, S).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a tricarboxylic acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules that contain atoms outside a set of typical organic elements.
    # Allowed elements: H, C, N, O, S.
    allowed_atomic_nums = {1, 6, 7, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non-organic atom: {atom.GetSymbol()} (atomic num: {atom.GetAtomicNum()})"

    # Define two SMARTS patterns for carboxyl groups:
    # One for protonated carboxylic acid -C(=O)[OH] (the -OH oxygen is sp2, has 2 neighbors, and one hydrogen)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    # One for deprotonated carboxylate group -C(=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    if acid_pattern is None or carboxylate_pattern is None:
        return False, "Error creating SMARTS patterns for carboxyl groups"
    
    # Get substructure matches for both patterns.
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # We use the index of the carbon atom in the carboxyl group (first atom in the match)
    # to count unique carboxyl groups.
    carboxyl_carbon_indices = set()
    for match in acid_matches:
        carboxyl_carbon_indices.add(match[0])
    for match in carboxylate_matches:
        carboxyl_carbon_indices.add(match[0])
    
    n_carboxyl_groups = len(carboxyl_carbon_indices)
    
    if n_carboxyl_groups != 3:
        return False, f"Found {n_carboxyl_groups} carboxyl group(s); expected exactly 3"
    
    # If we reach here, the molecule has exactly three carboxyl groups and is free of disallowed atoms.
    return True, "Contains exactly three carboxyl groups and no disallowed atoms, consistent with a tricarboxylic acid"

# Example usage (uncomment below lines for testing):
# test_smiles = [
#     "C[C@H](\\C=C\\C=C(\\C)[C@H]1CN[C@@H]([C@H]1CC(O)=O)C(O)=O)C(O)=O",  # domoic acid (true acid)
#     "O[C@@H]([C@H](CC(O)=O)C(O)=O)C(O)=O",  # D-erythro-isocitric acid (true acid)
#     "OC(=O)CCP(CCC(O)=O)CCC(O)=O",  # TCEP (contains P, should be excluded)
#     "CC(=O)N[C@@H](CC([O-])=O)C(=O)N[C@@H](CCC([O-])=O)C([O-])=O"  # Ac-Asp-Glu peptide fragment (false positive)
# ]
# for s in test_smiles:
#     result, reason = is_tricarboxylic_acid(s)
#     print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")