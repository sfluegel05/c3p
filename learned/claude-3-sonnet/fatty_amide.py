"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:33233 fatty amide
A monocarboxylic acid amide derived from a fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        return False, "Invalid SMILES string"
    
    # Look for amide pattern (C(=O)N[!#1])
    amide_pattern = Chem.MolFromSmarts("C(=O)N[!#1]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"
    
    # Look for linear carbon chain of at least 5 carbon atoms
    chain_pattern = Chem.MolFromSmarts("[CX4][CX4]~[CX4]~[CX4]~[CX4]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No linear carbon chain of at least 5 atoms found"
    
    # Check if the carbon chain is not part of a ring system
    ring_info = mol.GetRingInfo()
    for match in chain_matches:
        is_ring = any(ring_info.IsBondInRingOfSize(mol.GetBondBetweenAtoms(match[i], match[i+1]).GetIdx(), 3) for i in range(len(match)-1))
        if is_ring:
            return False, "Carbon chain is part of a ring system"
    
    # Count rotatable bonds and molecular weight
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    if n_rotatable < 5:
        return False, "Not enough rotatable bonds for a fatty amide"
    if mol_wt < 200:
        return False, "Molecular weight too low for a fatty amide"
    
    return True, "Molecule contains an amide group and a linear carbon chain of at least 5 atoms"