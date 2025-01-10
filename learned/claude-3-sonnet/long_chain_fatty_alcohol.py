"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
"""
Classifies: long-chain fatty alcohol
Definition: A fatty alcohol with a chain length ranging from C13 to C22
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"
    
    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 13:
        return False, f"Carbon chain too short (C{carbon_count}, need C13-C22)"
    if carbon_count > 22:
        return False, f"Carbon chain too long (C{carbon_count}, need C13-C22)"
        
    # Check if molecule is primarily a hydrocarbon chain
    # Count proportion of carbons and hydrogens vs other atoms
    total_atoms = mol.GetNumAtoms()
    h_count = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    ch_proportion = (carbon_count + h_count) / (total_atoms + h_count)
    
    if ch_proportion < 0.75:  # At least 75% should be C and H
        return False, "Not primarily a hydrocarbon structure"
    
    # Check for excessive rings - fatty alcohols typically have 0-1 rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count > 1:
        return False, f"Too many rings ({ring_count}) for a fatty alcohol"
        
    # Check for longest carbon chain
    # Convert to alkane by replacing all bonds with single bonds
    mol_alkane = Chem.RWMol(mol)
    for bond in mol_alkane.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
    
    # Find longest carbon chain
    longest_chain = 0
    for atom in mol_alkane.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            paths = Chem.FindAllPathsOfLengthN(mol_alkane, 20, useBonds=False, useHs=False)
            for path in paths:
                chain_length = sum(1 for idx in path if mol_alkane.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                longest_chain = max(longest_chain, chain_length)
    
    if longest_chain < 13:
        return False, f"Longest continuous carbon chain too short (C{longest_chain}, need C13-C22)"
        
    # Additional check for aromatic character
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    if aromatic_atoms > 6:  # Allow one benzene ring
        return False, "Too many aromatic atoms for a fatty alcohol"
    
    return True, f"Contains hydroxyl group with appropriate carbon chain length (C{carbon_count})"