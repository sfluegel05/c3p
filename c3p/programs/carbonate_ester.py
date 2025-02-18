"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester typically has the functional group O=C(O-)O- attached to organic groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to get molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a comprehensive set of SMARTS patterns for carbonate esters
    patterns = [
        Chem.MolFromSmarts("[$([C]=[O])]-[$([O]C)]"), # linear carbonate ester O=C-O-C
        Chem.MolFromSmarts("[$([C]-O-[C](=O)-[O]-C)]"), # cyclic or connected structures
        
        Chem.MolFromSmarts("O=C(OC)OC"), # typical dialkyl carbonate ester
        Chem.MolFromSmarts("O=C1OC1"),   # small cyclic carbonate (e.g., ethylene carbonate)
        Chem.MolFromSmarts("[$([O]C(=O)O)]")  # generic cyclic carbonate pattern
    ]

    # Check if the molecule matches any known carbonate ester pattern
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbonate ester functional group"
    
    return False, "No carbonate ester functional group found"