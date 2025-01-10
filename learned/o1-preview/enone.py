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
    
    # SMARTS pattern for alpha,beta-unsaturated ketone (enone)
    # This pattern matches a C=C double bond connected to a C=O group, where the carbonyl carbon has no hydrogens attached
    enone_pattern = Chem.MolFromSmarts('C=CC(=O)[C;H0]')
    
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains alpha,beta-unsaturated ketone (enone) group with R(4) ≠ H"
    else:
        return False, "No alpha,beta-unsaturated ketone (enone) group found"