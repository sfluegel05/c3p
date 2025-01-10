"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    Gangliosides are glycosphingolipids that contain one or more sialic acid residues
    linked on the sugar chain, and a ceramide backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined patterns for ceramide backbone
    # Ceramide is more varied; include the amide linkage and long chain portions
    ceramide_pattern = Chem.MolFromSmarts("NC(=O)C[C@H](O)CO")  # Including different chain lengths and positions
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No valid ceramide backbone found"

    # Refined patterns for sialic acid residues (using "C[C@H](O1)C(=O)O[C@@H]1" part more flexibly)
    # Generalized sialic acid indicator by focusing more on characteristic linkage
    sialic_acid_pattern = Chem.MolFromSmarts("[C@H](O[C@@H]([C@@H](CO)O)[C@H](NC(=O)C))")
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No valid sialic acid residue found"

    # Relaxed glycan pattern to efficiently match glucose and galactose configurations
    # Focusing on common glycosidic linkage patterns
    glycan_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@@H](CO)O[C@H]1")  # Flexible glucose/galactose pattern
    glycan_matches = mol.GetSubstructMatches(glycan_pattern)
    if len(glycan_matches) < 1:
        return False, "Missing necessary glycan residues"
        
    return True, "Contains ceramide backbone with sialic acid and necessary glycan residues"