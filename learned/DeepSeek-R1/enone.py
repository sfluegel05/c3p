"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:39176 enone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with the structure R1R2C=CR3-C(=O)R4 where R4 != H.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for enone: conjugated C=C adjacent to a ketone group
    enone_pattern = Chem.MolFromSmarts("[CX3]=[CX3]-[CX3](=[OX1])([#6])[#6]")
    
    # Check for matches
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains conjugated alpha,beta-unsaturated ketone"
    else:
        return False, "No conjugated alpha,beta-unsaturated ketone found"