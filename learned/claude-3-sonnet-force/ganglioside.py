"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: CHEBI:17619 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid with one or more sialic acids linked on the sugar chain.

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

    # Look for ceramide backbone (long alkyl chain + amide)
    ceramide_pattern = Chem.MolFromSmarts("CCCCC(=O)NC")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"
    
    # Look for N-acetylneuraminic acid (Neu5Ac) residue(s)
    neu5ac_pattern = Chem.MolFromSmarts("C[C@@]([C@@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)NC(C)=O)([C@@H]([C@@H](CO)O)O)(O)C(=O)O")
    neu5ac_matches = mol.GetSubstructMatches(neu5ac_pattern)
    if not neu5ac_matches:
        return False, "No N-acetylneuraminic acid (Neu5Ac) residues found"
    
    # Look for glycosidic linkages between sugar residues
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H][C@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if not glycosidic_matches:
        return False, "No glycosidic linkages found"
    
    # Count heavy atoms - gangliosides are typically large molecules
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 50:
        return False, "Too few heavy atoms for ganglioside"
    
    return True, "Contains ceramide backbone with N-acetylneuraminic acid (Neu5Ac) residue(s) and glycosidic linkages"