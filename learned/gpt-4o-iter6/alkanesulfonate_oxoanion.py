"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by a distinctive sulfonate group attached
    to an alkyl chain (R) represented by the pattern R-CS([O-])(=O)=O.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced sulfonate pattern to ensure connectivity context
    sulfonate_pattern = Chem.MolFromSmarts("C[S](=O)(=O)[O-]")
    
    # Check if the molecule contains the sulfonate group
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "Does not contain the characteristic alkanesulfonate sulfonate group"
    
    # Validate the presence of the alkane chain connected to the sulfonate
    # Consider adjacent carbons forming a linear or branched pattern
    # {contextual validation logic can go here if further specificity is needed}
    
    return True, "Contains the characteristic alkanesulfonate sulfonate group"

# Example use-case, add test cases to verify correctness and adjust if needed