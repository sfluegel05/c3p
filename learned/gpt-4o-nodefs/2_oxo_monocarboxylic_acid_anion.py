"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is classified as a 2-oxo monocarboxylic acid anion
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for strict carboxylate and 2-oxo group pattern (O=C-[C=O][O-])
    carboxylate_with_oxo_pattern = Chem.MolFromSmarts("O=C-[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_with_oxo_pattern):
        return False, "No 2-oxo and carboxylate linkage found"

    # Initiate a more specific search given structural complexity
    carb_position_oxo_pattern = Chem.MolFromSmarts("C(=O)[CX3](=O)C([O-])")
    if not mol.HasSubstructMatch(carb_position_oxo_pattern):
        return False, "No valid 2-oxo and monocarboxylate anion configuration"

    # Check for only one such specific contiguous pattern
    pattern_matches = mol.GetSubstructMatches(carboxylate_with_oxo_pattern)
    if len(pattern_matches) != 1:
        return False, f"Pattern mismatch; contains non-matching configurations"

    # Additional molecular checks to ensure it's the main anionic group
    total_neg_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    if total_neg_charge != -1:
        return False, "Total charge does not match expected single charge"

    return True, "Molecule matches 2-oxo monocarboxylic acid anion structure"

# Example usage
example_smiles = "NCCCC(=O)C([O-])=O"  # 5-amino-2-oxopentanoate
result, reason = is_2_oxo_monocarboxylic_acid_anion(example_smiles)
print(f"Classification result: {result}, Reason: {reason}")