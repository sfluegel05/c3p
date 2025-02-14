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
    
    # Look for nitrile (-C#N) and alcohol (-OH) groups
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#N")
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(nitrile_pattern) or not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Missing nitrile and/or alcohol groups"
    
    # Look for the cyanohydrin substructure: R-C(OH)(C#N)-R'
    cyanohydrin_pattern = Chem.MolFromSmarts("[CX4]([OH])([CX2]#N)")
    if not mol.HasSubstructMatch(cyanohydrin_pattern):
        return False, "Missing cyanohydrin substructure"
    
    # Check that the -OH and -C#N groups are attached to the same carbon
    cyanohydrin_atoms = mol.GetSubstructMatches(cyanohydrin_pattern)[0]
    carbon_atom = mol.GetAtomWithIdx(cyanohydrin_atoms[0])
    nitrile_atom = carbon_atom.GetNeighbors()[0]
    alcohol_atom = carbon_atom.GetNeighbors()[1]
    if nitrile_atom.GetAtomicNum() != 7 or alcohol_atom.GetAtomicNum() != 8:
        return False, "Nitrile and alcohol groups not attached to the same carbon"
    
    return True, "Contains the cyanohydrin substructure: R-C(OH)(C#N)-R'"