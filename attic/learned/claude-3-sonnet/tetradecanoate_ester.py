"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
"""
Classifies: tetradecanoate ester
A fatty acid ester obtained by condensation of the carboxy group of tetradecanoic acid 
with a hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule contains a tetradecanoate ester group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for non-cyclic ester group pattern
    # [C:1] must be sp3 carbon (no double bonds)
    # Must not be in ring
    ester_pattern = Chem.MolFromSmarts("[CH3X4:1][CH2X4:2][CH2X4:3][CH2X4:4][CH2X4:5][CH2X4:6][CH2X4:7][CH2X4:8][CH2X4:9][CH2X4:10][CH2X4:11][CH2X4:12][CH2X4:13][CX3:14](=[OX1:15])[OX2:16][*:17]")
    
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No tetradecanoate ester group found"
        
    matches = mol.GetSubstructMatches(ester_pattern)
    
    for match in matches:
        # Get the matched atoms
        chain_atoms = match[0:13]  # First 13 carbons
        carbonyl_carbon = match[13]
        ester_oxygen = match[15]
        
        # Verify none of the chain atoms are in a ring
        ring_atoms = mol.GetRingInfo().AtomRings()
        if any(atom_idx in {atom for ring in ring_atoms for atom in ring} for atom_idx in chain_atoms):
            continue
            
        # Verify chain carbons don't have additional carbon substituents
        has_branches = False
        for atom_idx in chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            carbon_neighbors = len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6])
            if carbon_neighbors > 2:  # More than 2 carbon neighbors means branching
                has_branches = True
                break
                
        if has_branches:
            continue
            
        # Verify no modifications on the chain (OH groups, double bonds, etc)
        has_modifications = False
        for atom_idx in chain_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetTotalNumHs() + len([n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]) < 4:
                has_modifications = True
                break
                
        if has_modifications:
            continue
            
        # If we get here, we've found a valid tetradecanoate ester group
        return True, "Contains tetradecanoate (myristoyl) ester group"
        
    return False, "No valid tetradecanoate ester group found"