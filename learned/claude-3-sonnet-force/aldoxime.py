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

    # Check if oxime is part of the main molecular scaffold
    scaffold_atoms = set(range(mol.GetNumAtoms()))
    for match in oxime_matches:
        atom_idx = match[0]
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = atom.GetNeighbors()
        scaffold_atoms.update(neighbor.GetIdx() for neighbor in neighbors)

    if len(scaffold_atoms) != mol.GetNumAtoms():
        return False, "Oxime group is not part of the main molecular scaffold"

    # Check for exceptional cases (e.g., oximes attached to heteroatoms or in conjugated systems)
    exceptional_pattern = Chem.MolFromSmarts("[CH]=[N][OH]~[cX3]")
    exceptional_matches = mol.GetSubstructMatches(exceptional_pattern)
    if exceptional_matches:
        return True, "Contains oxime functional group in an exceptional environment"

    # Check for common aldoxime environments
    aldoxime_patterns = [
        Chem.MolFromSmarts("[CH]=[N][OH]~[C](=O)"), # Acyl oximes
        Chem.MolFromSmarts("[CH]=[N][OH]~[C](-[OH])=O"), # Hydroxy aldoximes
        Chem.MolFromSmarts("[CH]=[N][OH]~[cX3](-[#6])~[#6]"), # Conjugated aldoximes
        Chem.MolFromSmarts("[CH]=[N][OH]~[cX3](-[#7])~[#7]"), # Heteroatom-containing aldoximes
        # Add more patterns as needed
    ]
    for pattern in aldoxime_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains aldoxime functional group in a common environment"

    return False, "Oxime group is not in a typical aldoxime environment"