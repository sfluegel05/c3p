"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid component esterified with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Butyric acid ester linkage pattern
    # The pattern allows for variations by considering potential isomers and flexibility in the ester component
    # Match a C3 or longer chain following "CCC(=O)O" that allows for alkyl or aryl variations following the ester group
    pat1 = Chem.MolFromSmarts("CCC(=O)OC")
    pat2 = Chem.MolFromSmarts("C[CH](C)C(=O)O")  # Including possible isomers like isobutyrate
    pat3 = Chem.MolFromSmarts("CC(C)C(=O)O")     # Including tert-butyrate

    # Check for match with any defined pattern
    if mol.HasSubstructMatch(pat1) or mol.HasSubstructMatch(pat2) or mol.HasSubstructMatch(pat3):
        return True, "Contains butyric acid ester linkage"

    return False, "No butyric acid ester linkage found"