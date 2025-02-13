"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    Gangliosides are glycosphingolipids that contain one or more sialic acid residues
    linked on the sugar chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for ceramide backbone pattern
    ceramide_pattern = Chem.MolFromSmarts("NC([C@@H](O)CO)C(=O)")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Look for sialic acid residue (NeuAc)
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@H](O1)C(=O)O[C@@H]([C@@H](CO)O)[C@@H]1NC=O")
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No sialic acid residue found"
        
    # Look for glycan residues: glucose and galactose
    glycan_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](CO)O[C@H]1O[C@H]1[C@H](O)[C@H](C=O)O[C@@H]1O")
    glycan_matches = mol.GetSubstructMatches(glycan_pattern)
    if len(glycan_matches) < 1:
        return False, "Missing necessary glycan residues"
        
    return True, "Contains ceramide backbone with sialic acid and glycan residues"