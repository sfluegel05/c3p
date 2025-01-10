"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone with a 1,3-diphenylpropenone structure (ArCH=CH(=O)Ar) and its derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible α,β-unsaturated ketone pattern
    alpha_beta_unsat_ketone_pattern = Chem.MolFromSmarts("C=CC(=O)[#6]")
    if not mol.HasSubstructMatch(alpha_beta_unsat_ketone_pattern):
        return False, "No α,β-unsaturated carbonyl group found"

    # Check specifically for the aryl-ketone-aryl structure (1,3-diphenylpropenone)
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)c2ccccc2")
    partial_chalcone_pattern = Chem.MolFromSmarts("cC=CC=O")
    
    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Contains the primary 1,3-diphenylpropenone structure"
    elif mol.HasSubstructMatch(partial_chalcone_pattern):
        return True, "Contains a variant structure similar to chalcone with modified aryl components"

    return False, "Does not contain a chalcone structure"