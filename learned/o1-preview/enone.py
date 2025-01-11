"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:xxxxx enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is defined as an alpha,beta-unsaturated ketone of general formula R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) ≠ H),
    in which the C=O function is conjugated to a C=C double bond at the alpha,beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for enones in both orientations
    enone_patterns = [
        Chem.MolFromSmarts('[C;H0](=O)-[C]=[C]'),  # O=C-C=C
        Chem.MolFromSmarts('[C]=[C]-[C;H0](=O)')   # C=C-C=O
    ]
    
    # Check if the molecule matches any enone pattern
    for pattern in enone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains alpha,beta-unsaturated ketone (enone) group with R(4) ≠ H"
    
    # No enone patterns matched
    return False, "No alpha,beta-unsaturated ketone (enone) group found"