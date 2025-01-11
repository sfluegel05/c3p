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

    # More flexible pattern to catch common features of 2-oxo monocarboxylic acid anions
    oxo_ketone_pattern = Chem.MolFromSmarts("[CX3]=O")  # C=O group
    carboxylate_anion_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # C(=O)[O-] group

    # Check for oxo (ketone) groups
    if not mol.HasSubstructMatch(oxo_ketone_pattern):
        return False, "No 2-oxo (ketone) group found"

    # Check for carboxylate anion groups
    if not mol.HasSubstructMatch(carboxylate_anion_pattern):
        return False, "No monocarboxylate anion group found"

    return True, "Molecule matches 2-oxo monocarboxylic acid anion structure"

# Example usage
example_smile = "NCCCC(=O)C([O-])=O"  # 5-amino-2-oxopentanoate
result, reason = is_2_oxo_monocarboxylic_acid_anion(example_smile)
print(f"Classification result: {result}, Reason: {reason}")