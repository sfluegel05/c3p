"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide contains a macrocyclic lactone ring with 12 or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify large rings (12 or more members)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    large_rings = [ring for ring in ring_info.AtomRings() if len(ring) >= 12]
    
    if not large_rings:
        return False, "No ring with 12 or more members found"

    # Check for a macrocyclic lactone (cyclic ester) presence
    # Improved SMARTS pattern for macrocyclic lactone:
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~1")
    
    for ring in large_rings:
        submol = Chem.PathToSubmol(mol, ring)
        if submol.HasSubstructMatch(lactone_pattern):
            return True, "Contains a macrocyclic lactone ring with 12 or more members"

    return False, "No macrocyclic lactone structure found in large rings"