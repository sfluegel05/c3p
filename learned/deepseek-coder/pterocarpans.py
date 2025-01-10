"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: CHEBI:XXXXX pterocarpans
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    A pterocarpan has a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core pterocarpan skeleton pattern
    # This pattern captures the fused benzofuro[3,2-c]chromene system
    pterocarpan_pattern = Chem.MolFromSmarts("[O;R]1[C;R]2[C;R]([C;R]=3[C;R]1=[C;R][C;R]=[C;R]3)[C;R][O;R][C;R]4=[C;R]2[C;R]=[C;R][C;R]=[C;R]4")
    
    # Check if the molecule matches the pterocarpan skeleton
    if not mol.HasSubstructMatch(pterocarpan_pattern):
        return False, "Molecule does not match the pterocarpan skeleton"

    # Additional checks for typical pterocarpan features
    # Count the number of aromatic rings (should be at least 2)
    aromatic_rings = Chem.GetSSSR(mol)
    if len(aromatic_rings) < 2:
        return False, "Molecule does not have enough aromatic rings"

    # Check for the presence of oxygen atoms (at least 2)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Molecule does not have enough oxygen atoms"

    # Check for the presence of a benzofuran and chromene system
    benzofuran_pattern = Chem.MolFromSmarts("[O;R]1[C;R]2=[C;R][C;R]=[C;R][C;R]=[C;R]2[C;R]1")
    chromene_pattern = Chem.MolFromSmarts("[O;R]1[C;R]2=[C;R][C;R]=[C;R][C;R]=[C;R]2[C;R]1")
    
    if not mol.HasSubstructMatch(benzofuran_pattern) or not mol.HasSubstructMatch(chromene_pattern):
        return False, "Molecule does not contain both benzofuran and chromene systems"

    return True, "Molecule matches the pterocarpan skeleton with benzofuran and chromene systems"