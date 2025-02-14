"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    Azoles are monocyclic 5-membered heterocycles containing at least one nitrogen and
    possibly other heteroatoms (N, S, O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find 5-membered rings with at least 1 Nitrogen
    azole_pattern = Chem.MolFromSmarts("[n;R5][n,o,s;R5][*;R5][*;R5][*;R5]")
    
    ring_matches = mol.GetSubstructMatches(azole_pattern)

    if not ring_matches:
          return False, "No azole ring found"
    
    # Check the ring is monocyclic (exactly 1 ring)
    if Chem.GetSSSR(mol) != 1:
       return False, "Found fused ring or no ring"


    return True, "Contains a monocyclic 5-membered heterocycle with at least one nitrogen."