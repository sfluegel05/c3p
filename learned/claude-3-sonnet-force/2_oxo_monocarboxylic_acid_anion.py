"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:24351 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    A 2-oxo monocarboxylic acid anion has an oxo group at the 2-position and a negatively charged carboxylate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylate group [-C(=O)[O-]]
    carboxylate_pattern = Chem.MolFromSmarts("[C-](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Found {len(carboxylate_matches)} carboxylate groups, need exactly 1"
    
    # Look for oxo group [C(=O)-]
    oxo_pattern = Chem.MolFromSmarts("[C](=O)(-*)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) != 1:
        return False, f"Found {len(oxo_matches)} oxo groups, need exactly 1"
    
    # Check if oxo group is at 2-position relative to carboxylate
    carboxylate_atom = mol.GetAtomWithIdx(carboxylate_matches[0][0])
    oxo_atom = mol.GetAtomWithIdx(oxo_matches[0][0])
    
    # Find the shortest path between the carboxylate carbon and the oxo carbon
    path = Chem.rdmolops.GetShortestPath(mol, carboxylate_atom.GetIdx(), oxo_atom.GetIdx())
    
    # Check if the path length is 2 (i.e., oxo group at 2-position)
    if len(path) != 3:
        return False, "Oxo group not at 2-position relative to carboxylate"
    
    # Check if the bond between the carboxylate carbon and the adjacent atom is a single bond
    if mol.GetBondBetweenAtoms(carboxylate_atom.GetIdx(), path[1]).GetBondType() != BondType.SINGLE:
        return False, "Carboxylate group not attached via a single bond"
    
    # Check for only one carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[C](=O)[O;!-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) > 0:
        return False, "Contains non-ionized carboxylic acid groups"
    
    return True, "Contains a carboxylate group and an oxo group at the 2-position"