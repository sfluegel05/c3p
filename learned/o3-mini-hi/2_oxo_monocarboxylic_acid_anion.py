"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-oxo monocarboxylic acid anion
Definition: An oxo monocarboxylic acid anion in which the oxo (carbonyl) group is located at the 2‑position.
This generally corresponds to α‑keto acid anions with a core structure R–C(=O)–C(=O)[O–].
"""

from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    That is, the molecule must contain an α-keto acid motif where the carboxylate group [-C(=O)[O-]]
    is directly attached to a carbonyl (oxo) group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for the classification decision
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the 2-oxo monocarboxylic acid anion motif.
    # The pattern "[#6X3](=O)[#6X3](=O)[O-]" matches:
    #   - a sp2-hybridized carbon (#6X3) with a carbonyl (C=O)
    #   - directly bonded to another sp2-hybridized carbon that forms a carbonyl (C=O)
    #     and is bound to an oxygen with a negative charge ([O-]), i.e. the carboxylate.
    motif = Chem.MolFromSmarts("[#6X3](=O)[#6X3](=O)[O-]")
    
    if mol.HasSubstructMatch(motif):
        return True, "Contains the 2-oxo monocarboxylic acid anion (α-keto acid anion) motif."
    else:
        return False, "Does not contain the 2-oxo monocarboxylic acid anion motif."
        
# Examples (uncomment the lines below to test the function):
# test_smiles = "CC(=O)C([O-])=O"  # pyruvate, a classic example of an α-keto acid anion
# result, reason = is_2_oxo_monocarboxylic_acid_anion(test_smiles)
# print(result, reason)