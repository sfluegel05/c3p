"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: phenyl acetates (CHEBI:xxxxx)
An acetate ester obtained by formal condensation of the carboxy group of acetic acid with the hydroxy group of any phenol.
"""
from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate has an acetate ester group (O-C(=O)-CH3) attached to a benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenyl acetate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for acetate ester attached to a benzene ring
    # [c&r6] is an aromatic carbon in a 6-membered ring (benzene)
    # connected to O-C(=O)-CH3
    pattern = Chem.MolFromSmarts("[c&r6]-OC(=O)C")
    
    # Check for matches
    if not mol.HasSubstructMatch(pattern):
        return False, "No acetate ester group attached to a benzene ring"
    
    return True, "Contains acetate ester group attached to a benzene ring"