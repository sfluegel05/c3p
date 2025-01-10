"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is defined as 4-(2-aminoethyl)pyrocatechol and derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catecholamine core pattern
    # Benzene ring with hydroxyls at positions 1 and 2, aminoethyl chain at position 4
    catecholamine_pattern = Chem.MolFromSmarts("""
    [#6]-1
      -[#6](=[#6]-[#6]-[#6]-[#6]-1-)
      -[#8]-[#1]
      -[#6]-[#6]-[#7]
    """)
    
    # Alternative SMARTS pattern with explicit positions
    # c1c(O)cc(O)c(c1)CCN
    catecholamine_pattern = Chem.MolFromSmarts('c1c(O)cc(O)c(c1)CCN')
    
    if not mol.HasSubstructMatch(catecholamine_pattern):
        return False, "Molecule does not match catecholamine core structure"
    
    return True, "Molecule matches catecholamine core structure"