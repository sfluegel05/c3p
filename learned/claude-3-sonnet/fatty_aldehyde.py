"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde has an aldehyde group at one end of a carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules with rings (except small rings that might be errors in SMILES)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 0:
        # Check ring sizes - only allow 3-membered rings which might be SMILES artifacts
        ssr = Chem.GetSymmSSSR(mol)
        for ring in ssr:
            if len(ring) > 3:
                return False, "Contains ring structure - not a fatty aldehyde"

    # Look for aldehyde groups (-CH=O)
    aldehyde_pattern = Chem.MolFromSmarts("[CH1,CH0]([H,#1,*])=O")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    if not aldehyde_matches:
        return False, "No aldehyde group found"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4:
        return False, "Carbon chain too short for fatty aldehyde"

    # Check for carboxylic acids
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains carboxylic acid group"

    # Check for linear carbon chain
    chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficient carbon chain found"

    # Check that most atoms are carbons (allowing for some heteroatoms)
    total_heavy_atoms = mol.GetNumHeavyAtoms()
    if carbon_count < total_heavy_atoms * 0.7:
        return False, "Too many heteroatoms for fatty aldehyde"

    # Check for excessive branching
    branching_points = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) > 2:
                branching_points += 1
    
    if branching_points > 2:  # Allow some branching but not too much
        return False, "Too many branch points for fatty aldehyde"

    # Success - determine if saturated or unsaturated
    double_bond_pattern = Chem.MolFromSmarts("[C]=[C]")
    has_double_bonds = mol.HasSubstructMatch(double_bond_pattern)
    
    if has_double_bonds:
        return True, "Unsaturated fatty aldehyde with terminal aldehyde group"
    else:
        return True, "Saturated fatty aldehyde with terminal aldehyde group"