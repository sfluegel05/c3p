"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester
Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
This means that one of the ester groups in the molecule has an acyl part corresponding to n-butyrate (CH3CH2CH2COO-),
typically represented in SMILES as the fragment "OC(=O)CCC".
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    The function checks for the presence of an ester moiety where the acyl component matches that of butyric acid.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a butyrate ester fragment, False otherwise.
        str: A reason for the classification result.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern that represents the butyrate ester fragment.
    # "O[C](=O)CCC" matches an oxygen (ester linkage) attached to a carbonyl carbon (C(=O))
    # that is bonded to three aliphatic carbons (CCC), corresponding to n-butyrate.
    butyrate_smarts = "O[C](=O)CCC"
    butyrate_pattern = Chem.MolFromSmarts(butyrate_smarts)
    
    # Check if the molecule contains the butyrate ester fragment.
    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Found butyrate ester fragment (O-C(=O)CCC) in the molecule."
    else:
        return False, "No butyrate ester fragment (O-C(=O)CCC) found in the molecule."