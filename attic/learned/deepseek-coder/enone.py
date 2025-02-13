"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with the general formula R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) =/= H).

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

    # Define the enone pattern: C=C-C=O where the carbonyl carbon is not bonded to hydrogen
    # More flexible pattern that allows for any number of atoms between the double bond and carbonyl
    enone_pattern = Chem.MolFromSmarts("[CX3]=[CX3][CX3]=[OX1]")
    
    # Check for the pattern
    if mol.HasSubstructMatch(enone_pattern):
        return True, "Contains alpha,beta-unsaturated ketone (enone) structure"
    
    # Alternative pattern for cases where there might be conjugation through aromatic systems
    aromatic_enone_pattern = Chem.MolFromSmarts("[c]:[c][CX3]=[OX1]")
    if mol.HasSubstructMatch(aromatic_enone_pattern):
        return True, "Contains aromatic alpha,beta-unsaturated ketone (enone) structure"

    return False, "No alpha,beta-unsaturated ketone pattern found"