"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    Gangliosides are glycosphingolipids that contain one or more sialic acid residues
    linked to the sugar chain, and a ceramide backbone.

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

    # Generalized pattern for ceramide backbone
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@H](CO)C([O])[C@H][C@@H]")  # Core amide and part of sphingosine structure
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No valid ceramide backbone found"

    # Generalized pattern for sialic acid residues
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@H](C(=O)[O-])O[C@H]")  # Basic sialic acid recognition with carboxylate
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No valid sialic acid residue found"

    # Checking for a minimum number of sugar residues (assuming glycan includes glucose/galactose)
    glycan_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")  # Pattern for common glucose/galactose units
    glycan_matches = mol.GetSubstructMatches(glycan_pattern)
    if len(glycan_matches) < 3:
        return False, "Missing necessary glycan residues"
        
    return True, "Contains ceramide backbone with sialic acid and necessary glycan residues"