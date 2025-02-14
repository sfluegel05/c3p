"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: CHEBI:51856 tricarboxylic acid
"An oxoacid containing three carboxy groups"
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Remove counterions and neutralize charges
    mol = AllChem.RemoveHs(mol, implicitOnly=True, updateExplicitCount=True)
    AllChem.Uncharger(mol)
    
    # Look for 3 carboxy groups (-C(=O)O-)
    carboxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1]")
    carboxy_matches = mol.GetSubstructMatches(carboxy_pattern)
    if len(carboxy_matches) != 3:
        return False, f"Found {len(carboxy_matches)} carboxy groups, need exactly 3"
    
    # Look for oxo group (C=O)
    oxo_pattern = Chem.MolFromSmarts("C=O")
    oxo_match = mol.HasSubstructMatch(oxo_pattern)
    if not oxo_match:
        return False, "No oxo group found"
    
    # Check connectivity of carboxy groups
    carboxy_atoms = set([atom.GetIdx() for match in carboxy_matches for atom in match])
    connected_carboxy_atoms = set()
    for atom_idx in carboxy_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        connected_atoms = [mol.GetAtomWithIdx(neighbor_idx).GetIdx() for neighbor_idx in atom.GetNeighbors()]
        connected_carboxy_atoms.update(connected_atoms)
    
    if len(connected_carboxy_atoms & carboxy_atoms) < 3:
        return False, "Carboxy groups not connected to the same core"
    
    # Check molecular weight
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 500:
        return False, "Molecular weight outside typical range for tricarboxylic acids"
    
    return True, "Molecule contains three carboxy groups and an oxo group, with carboxy groups connected to the same core"