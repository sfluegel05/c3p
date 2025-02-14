"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: CHEBI:35620 cyanohydrin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the formal addition of hydrogen cyanide to the C=O bond of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carbonyl group (C=O)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"
    
    # Look for cyanohydrin substructure: R-C(OH)(C#N)-R'
    cyanohydrin_pattern = Chem.MolFromSmarts("[CX4]([OH])([CX2]#N)")
    if not mol.HasSubstructMatch(cyanohydrin_pattern):
        return False, "Missing cyanohydrin substructure"
    
    # Check if the cyanohydrin substructure is formed by addition to a carbonyl
    cyanohydrin_atoms = mol.GetSubstructMatches(cyanohydrin_pattern)[0]
    carbon_atom = mol.GetAtomWithIdx(cyanohydrin_atoms[0])
    
    # Find the carbonyl carbon atom
    carbonyl_atoms = mol.GetSubstructMatches(carbonyl_pattern)[0]
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_atoms[0])
    
    # Check if the cyanohydrin carbon is the same as the carbonyl carbon
    if carbon_atom.GetIdx() != carbonyl_carbon.GetIdx():
        return False, "Cyanohydrin substructure not formed from addition to a carbonyl"
    
    return True, "Contains the cyanohydrin substructure formed by addition of hydrogen cyanide to a carbonyl group"