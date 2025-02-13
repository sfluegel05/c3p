"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester results from the esterification of the carboxy group of decanoic acid
    (which is a 10 carbon acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for esters specifically linked to a 10-carbon chain from decanoic acid
    ester_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)O")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Contains a decanoate ester linkage"
    
    return False, "No appropriate 10-carbon chain detected (decanoic acid derivative)"