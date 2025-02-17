"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: 3-hydroxy fatty acid:
Any fatty acid with a hydroxy functional group in the beta- or 3-position.
This function verifies that the molecule contains a free carboxylic acid group 
and that there is a beta-hydroxy motif (O=C(O)[CH2][CH](O)).
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid should contain a free carboxylic acid group, and a hydroxy
    (OH) substitution on the beta (3-) carbon relative to the carboxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, check for the presence of a carboxylic acid group.
    # The SMARTS "C(=O)[O;H]" looks for a carbon double-bonded to oxygen and attached to an OH.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid (COOH) group found, so not a fatty acid"
    
    # Next, check for the beta-hydroxy (3-hydroxy) motif.
    # The SMARTS "O=C(O)[CH2][CH](O)" represents a carboxylic acid (O=C(O))
    # followed by a CH2 (alpha carbon) and then a CH with an OH (beta carbon).
    beta_hydroxy_pattern = Chem.MolFromSmarts("O=C(O)[CH2][CH](O)")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        return False, "No beta-hydroxy group found relative to the carboxyl group"
    
    return True, "Contains a carboxylic acid group with a beta-hydroxy (3-hydroxy) substituent"
    
# Example test (uncomment to run if desired)
# test_smiles = "OCCCCCCCCCCCCCC[C@@H](O)CC(O)=O"  # (3R)-3,17-dihydroxymargaric acid
# result, reason = is_3_hydroxy_fatty_acid(test_smiles)
# print(result, reason)