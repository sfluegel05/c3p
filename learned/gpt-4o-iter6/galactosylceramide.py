"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a cerebroside with a galactose monosaccharide head group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a flexible sphingosine-like backbone pattern allowing variations
    sphingosine_pattern = Chem.MolFromSmarts("C[C@@H](O)C(N[C@H]1COC(O)C(O)C(O)C1O)")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No flexible sphingosine-like backbone found"
    
    # Look for amide linkage with long carbon chain
    amide_chain_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](C)C")
    if not mol.HasSubstructMatch(amide_chain_pattern):
        return False, "No amide linkage with long carbon chain found"
    
    # Look for galactose head group pattern - allow for variability like sulfation
    galactose_pattern = Chem.MolFromSmarts("C1(O[C@H](CO)[C@H](OC)C(O)C(O)C1O)")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No variable galactose head group found"

    return True, "Contains a flexible sphingosine-like backbone with amide linkage and galactose head group"