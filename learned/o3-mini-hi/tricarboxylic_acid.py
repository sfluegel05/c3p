"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: Tricarboxylic Acid 
Definition: An oxoacid containing three carboxy groups.
This version improves on the previous approach by:
  - Rejecting molecules that are not neutral (formal charge != 0)
  - Excluding molecules containing amide bonds (C(=O)N), since many false
    positives arose from peptides or protein fragments.
Then it detects exactly three unique carboxyl groups (in protonated or 
deprotonated form) using two SMARTS patterns.
Note: This classification is heuristic.
"""
from rdkit import Chem

def is_tricarboxylic_acid(smiles: str):
    """
    Determines whether a molecule is a tricarboxylic acid based on the SMILES string.
    
    A tricarboxylic acid is defined here as a neutral (formal charge zero) organic oxoacid
    containing exactly three carboxyl groups (–C(=O)[OH] or –C(=O)[O-]), and not having 
    amide bonds (which would suggest peptide/protein fragments).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tricarboxylic acid; False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES to an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules that contain atoms outside a set of allowed organic elements.
    # Allowed elements: H, C, N, O, S.
    allowed_atomic_nums = {1, 6, 7, 8, 16}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains non-organic atom: {atom.GetSymbol()} (atomic num: {atom.GetAtomicNum()})"
    
    # Check that molecule is neutral (formal charge must equal zero).
    if Chem.GetFormalCharge(mol) != 0:
        return False, "Molecule has nonzero formal charge; expected a neutral tricarboxylic acid"
    
    # Reject molecules that contain amide bonds.
    # The pattern "C(=O)N" is typical of an amide (peptide bond); its presence
    # suggests the molecule is a peptide or contains protein fragments.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if amide_pattern is not None and mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide bonds, indicative of peptides or protein fragments"
    
    # Define SMARTS patterns for the carboxyl (acid) groups.
    # One pattern for protonated carboxylic acid: -C(=O)[OH]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    # And one for deprotonated carboxylate: -C(=O)[O-]
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if acid_pattern is None or carboxylate_pattern is None:
        return False, "Error creating SMARTS patterns for carboxyl groups"
    
    # Get all substructure matches for both patterns.
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Count unique carboxyl groups by collecting the index of the carbon atom (first atom in the match).
    carboxyl_carbon_indices = set()
    for match in acid_matches:
        carboxyl_carbon_indices.add(match[0])
    for match in carboxylate_matches:
        carboxyl_carbon_indices.add(match[0])
    
    n_carboxyl_groups = len(carboxyl_carbon_indices)
    
    if n_carboxyl_groups != 3:
        return False, f"Found {n_carboxyl_groups} carboxyl group(s); expected exactly 3"
    
    return True, "Contains exactly three carboxyl groups, is neutral and free of amide bonds, consistent with a tricarboxylic acid"

# Example usage (uncomment for local testing):
# test_smiles = [
#     # True examples:
#     "C[C@H](\\C=C\\C=C(\\C)[C@H]1CN[C@@H]([C@H]1CC(O)=O)C(O)=O",  # domoic acid
#     "O[C@@H]([C@H](CC(O)=O)C(O)=O)C(O)=O",  # D-erythro-isocitric acid
#     # False examples:
#     "OC(=O)CCP(CCC(O)=O)CCC(O)=O",  # TCEP (has a P atom)
#     "CC(=O)N[C@@H](CC([O-])=O)C(=O)N[C@@H](CCC([O-])=O)C([O-])=O"  # peptide fragment with amide bonds
# ]
# for s in test_smiles:
#     result, reason = is_tricarboxylic_acid(s)
#     print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")