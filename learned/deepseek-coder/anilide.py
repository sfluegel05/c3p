"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: CHEBI:132109 anilide
Anilide is defined as any aromatic amide obtained by acylation of aniline.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide must contain an aromatic amide group where the nitrogen is directly attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the anilide pattern: aromatic ring with amide group attached
    # Pattern: [NX3]([CX3](=[OX1]))[c]
    anilide_pattern = Chem.MolFromSmarts("[NX3]([CX3](=[OX1]))[c]")
    
    # Check if the pattern exists in the molecule
    if not mol.HasSubstructMatch(anilide_pattern):
        return False, "No aromatic amide group found"
    
    # Get all matches of the pattern
    matches = mol.GetSubstructMatches(anilide_pattern)
    
    # Check if any match has nitrogen directly connected to aromatic carbon
    for match in matches:
        nitrogen_idx = match[0]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Check if nitrogen is connected to aromatic carbon
        for neighbor in nitrogen_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                return True, "Contains aromatic amide group with nitrogen attached to aromatic ring"
    
    return False, "No aromatic amide group with nitrogen attached to aromatic ring found"