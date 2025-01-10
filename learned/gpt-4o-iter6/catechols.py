"""
Classifies: CHEBI:33566 catechols
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol contains an o-diphenol component, which is a benzene ring with two adjacent hydroxyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a catechol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define catechol SMARTS pattern
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)cc(O)c1")
    
    # Check for catechol substructure
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains a catechol moiety (o-diphenol component)"
    else:
        return False, "No catechol moiety found"

# Examples and usage
examples = [
    "O[C@H]([C@H](OC(=O)\\C=C\\c1ccc(O)c(O)c1)C(O)=O)C(O)=O",
    "Oc1cc(O)cc(O)c1"
]
for example in examples:
    result, reason = is_catechols(example)
    print(f"SMILES: {example} -> Is Catechol: {result}, Reason: {reason}")