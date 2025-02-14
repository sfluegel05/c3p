"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:51856 tricarboxylic acid
"An oxoacid containing three carboxy groups"
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid is an oxoacid containing three carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove explicit hydrogens
    mol = Chem.RemoveHs(mol)
    
    # Look for 3 carboxy groups (-C(=O)O-)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if len(carboxy_matches) != 3:
        return False, f"Found {len(carboxy_matches)} carboxy groups, need exactly 3"
    
    # Check connectivity of carboxy groups
    connected_atoms = set()
    for match in carboxy_matches:
        for atom_idx in match:
            connected_atoms.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            connected_atoms.update(neighbor.GetIdx() for neighbor in atom.GetNeighbors())
    
    if len(connected_atoms) < 4:
        return False, "Carboxy groups not connected to the same core"
    
    # Check for oxo group (C=O)
    oxo_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No oxo group found"
    
    # Check molecular weight (optional)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 500:
        return False, "Molecular weight outside typical range for tricarboxylic acids"
    
    return True, "Molecule contains three carboxy groups and an oxo group, with carboxy groups connected to the same core"