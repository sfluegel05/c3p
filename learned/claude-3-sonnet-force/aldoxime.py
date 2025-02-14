"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:35781 aldoxime

An aldoxime is an oxime of an aldehyde, with the general structure RCH=NOH.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    An aldoxime is an oxime of an aldehyde, with the general structure RCH=NOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for oxime functional group (-CH=N-O-)
    oxime_pattern = Chem.MolFromSmarts("[CH]=[N][OH]")
    oxime_matches = mol.GetSubstructMatches(oxime_pattern)
    if not oxime_matches:
        return False, "No oxime functional group (-CH=N-O-) found"

    # Check if oxime is attached to an aldehyde carbon
    for match in oxime_matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetDegree() == 1 and atom.GetHybridization() == Chem.HybridizationType.SP2:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetDegree() == 1 and neighbor.GetAtomicNum() == 8:  # Check for C=O
                return True, "Contains aldoxime functional group (-CH=N-O-)"

    # Check for exceptional cases (e.g., oximes attached to heteroatoms or in conjugated systems)
    exceptional_pattern = Chem.MolFromSmarts("[CH]=[N][OH]~[cX3]")
    exceptional_matches = mol.GetSubstructMatches(exceptional_pattern)
    if exceptional_matches:
        return True, "Contains oxime functional group in an exceptional environment"

    return False, "Oxime not attached to aldehyde carbon or in an exceptional environment"