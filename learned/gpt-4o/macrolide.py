"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide has a macrocyclic lactone ring with 12 or more members.

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
    large_ring_exists = False

    for ring in ring_info.AtomRings():
        if len(ring) >= 12:
            large_ring_exists = True
            # Examine this specific ring for lactone characteristics
            submol = Chem.PathToSubmol(mol, ring)

            # Improved lactone search, allowing for heteroatoms and ensuring a ring is formed
            lactone_pattern = Chem.MolFromSmarts("[O;R][C;!H0;R](=O)")
            
            if submol.HasSubstructMatch(lactone_pattern):
                return True, "Contains a macrocyclic lactone ring with 12 or more members"

    if not large_ring_exists:
        return False, "No ring with 12 or more members found"

    return False, "No macrocyclic lactone structure found in large rings"