"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam compounds
A lactam in which the amide bond is contained within a four-membered ring, 
which includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Convert to neutral form if possible and add Hs
    mol = Chem.AddHs(mol)
    
    # More general SMARTS patterns for beta-lactam ring
    # Pattern 1: Basic beta-lactam core
    pattern1 = Chem.MolFromSmarts("[#7]1[#6][#6](=[#8])[#6]1")
    # Pattern 2: Alternative representation accounting for different bond types
    pattern2 = Chem.MolFromSmarts("[#7]1-[#6]-[#6](=O)-[#6]-1")
    # Pattern 3: Account for fused systems
    pattern3 = Chem.MolFromSmarts("[#7]1[#6][#6](=O)[#6]1[#6,#7,#8,#16]")
    
    patterns = [
        (pattern1, "basic beta-lactam core"),
        (pattern2, "alternative beta-lactam core"),
        (pattern3, "fused beta-lactam system")
    ]
    
    matches = []
    match_types = set()
    
    for pattern, desc in patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matches.extend(mol.GetSubstructMatches(pattern))
            match_types.add(desc)
    
    if not matches:
        return False, "No beta-lactam ring found"
    
    # Verify ring properties for each match
    ring_info = mol.GetRingInfo()
    valid_rings = 0
    
    for match in matches:
        ring_atoms = set(match)
        
        # Check if these atoms form a 4-membered ring
        for ring in ring_info.AtomRings():
            if len(set(ring).intersection(ring_atoms)) >= 4:
                # Verify presence of amide
                for atom_idx in ring:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomicNum() == 7:  # Nitrogen
                        # Check if nitrogen is connected to carbonyl
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 6:  # Carbon
                                for bond in neighbor.GetBonds():
                                    other_atom = bond.GetOtherAtom(neighbor)
                                    if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                                        valid_rings += 1
                                        break
    
    if valid_rings == 0:
        return False, "Found potential ring but missing required amide bond"
    
    # Additional validation for specific beta-lactam types
    found_types = ", ".join(match_types)
    return True, f"Contains {valid_rings} beta-lactam ring(s) ({found_types})"