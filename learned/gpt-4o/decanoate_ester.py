"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester results from the esterification of decanoic acid with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the decanoic acid part: a C10 chain with a terminal carboxylate
    decanoic_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)O")
    if not mol.HasSubstructMatch(decanoic_acid_pattern):
        return False, "No decanoic acid moiety found"

    # Look for the ester linkage: C(=O)O
    ester_linkage_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage found"

    return True, "Contains decanoic acid moiety with ester linkage"