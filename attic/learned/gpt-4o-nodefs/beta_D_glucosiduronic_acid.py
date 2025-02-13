"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule contains a beta-D-glucuronic acid moiety based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains beta-D-glucuronic acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string and check validity
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for beta-D-glucuronic acid moiety
    # 6-membered ring with OH groups and one COOH
    glucuronic_acid_smarts = "[C@H]1(O[C@H](CO)[C@@H]([C@H](O)[C@H]1O)C(=O)O)"
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)

    if mol.HasSubstructMatch(glucuronic_acid_pattern):
        return True, "Contains beta-D-glucuronic acid moiety"
    else:
        return False, "No beta-D-glucuronic acid moiety found"

# Test with examples
example_smiles = "O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)OC=2C=CC(=CC2)[N+]([O-])=O"
print(is_beta_D_glucosiduronic_acid(example_smiles))