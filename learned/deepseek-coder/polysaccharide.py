"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule consisting of large numbers of monosaccharide residues linked glycosidically.
    This term is commonly used only for those containing more than ten monosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a pattern for a monosaccharide unit (e.g., a hexose)
    monosaccharide_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@H]([C@@H]([C@H](O1)O)O)O)O)")
    
    # Find all matches of the monosaccharide pattern
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    
    # Count the number of monosaccharide units
    num_monosaccharides = len(monosaccharide_matches)
    
    # Check if the number of monosaccharide units is greater than 10
    if num_monosaccharides > 10:
        return True, f"Contains {num_monosaccharides} monosaccharide units, which is more than 10"
    else:
        return False, f"Contains {num_monosaccharides} monosaccharide units, which is not more than 10"