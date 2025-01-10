"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone moiety with at least one hydroxy group substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns - improved pattern for naphthoquinone framework
    naphthoquinone_pattern = Chem.MolFromSmarts('[#6]1[#6][#6]=[#6][#6][#6]2[#6](=O)[#6](=O)[#6][#6]12')
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    
    # Check for naphthoquinone structure (in various orientations)
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No naphthoquinone core found"
    
    # Check for at least one hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No hydroxy group found"
    
    # Ensure hydroxy group is part of naphthoquinone ring system
    for match in mol.GetSubstructMatches(naphthoquinone_pattern):
        naphthoquinone_atoms = set(match)
        hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
        
        # Check if any hydroxy is connected directly to the aromatic carbons in the naphthoquinone
        connected = False
        for h in hydroxy_matches:
            hydroxy_bonded_atoms = [bonded_atom.GetIdx() for bonded_atom in mol.GetAtomWithIdx(h[0]).GetNeighbors()]
            if any(atom_idx in naphthoquinone_atoms for atom_idx in hydroxy_bonded_atoms):
                connected = True
                break
        
        if connected:
            return True, "Contains naphthoquinone core with hydroxy group substitution"
    
    return False, "Hydroxy group not attached correctly to naphthoquinone"