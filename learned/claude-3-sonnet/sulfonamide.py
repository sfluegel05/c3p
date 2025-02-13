"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:22667 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.BRICS import BRICSDecompose

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is defined as an amide of a sulfonic acid with the general structure RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonamide substructure
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    if not mol.HasSubstructMatch(sulfonamide_pattern):
        return False, "No sulfonamide substructure found"
    
    # Check for carbon substituents on S and N
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur
            if not any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
                return False, "Sulfur not substituted with carbon"
        elif atom.GetAtomicNum() == 7:  # Nitrogen
            if not any(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
                return False, "Nitrogen not substituted with carbon"
    
    # Check for specific sulfonamide patterns
    brics_fragments = BRICSDecompose(mol)
    sulfonamide_patterns = ["c1ccccc1S(=O)(=O)N", "c1ccncc1S(=O)(=O)N", "C1CCCCC1S(=O)(=O)N"]
    for frag in brics_fragments:
        smarts = Chem.MolToSmarts(frag)
        if any(pat in smarts for pat in sulfonamide_patterns):
            return True, "Matched specific sulfonamide pattern"
    
    return True, "Contains sulfonamide substructure with carbon substituents"