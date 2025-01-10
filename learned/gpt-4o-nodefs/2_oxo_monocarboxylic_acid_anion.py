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

    # Define refined SMARTS patterns
    oxo_pattern = Chem.MolFromSmarts("C(=O)[CX4,CX3]")  # C group next to oxo
    monocarboxylate_anion_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # Carboxylate anion

    # Check for 2-oxo and carboxylate groups in proximal arrangement
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    anion_matches = mol.GetSubstructMatches(monocarboxylate_anion_pattern)
    
    # Ensure there is at least one oxo and one carboxylate and they are within 2 bonds
    for oxo in oxo_matches:
        for anion in anion_matches:
            if abs(oxo[0] - anion[0]) == 1:  # Ensure groups are adjacent
                return True, "Molecule matches 2-oxo monocarboxylic acid anion structure"
    
    return False, "Molecule does not match the structure of a 2-oxo monocarboxylic acid anion"

# Example usage
example_smile = "NCCCC(=O)C([O-])=O"  # 5-amino-2-oxopentanoate
result, reason = is_2_oxo_monocarboxylic_acid_anion(example_smile)
print(f"Classification result: {result}, Reason: {reason}")