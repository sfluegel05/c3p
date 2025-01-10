"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:35879 lactol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    Lactols are cyclic hemiacetals formed by intramolecular addition of a hydroxy 
    group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic cyclic hemiacetal pattern with sp3 carbon
    # [O;R] - ring oxygen
    # [C;R;X4] - sp3 carbon in ring with 4 connections
    # [OH1,OR] - hydroxy or alkoxy group
    # [!$(C=O)] - not part of a carbonyl
    hemiacetal_pattern = Chem.MolFromSmarts("[O;R][C;R;X4;!$(C=O)]([OH1,OR])")
    
    if not mol.HasSubstructMatch(hemiacetal_pattern):
        return False, "No lactol structure found"

    # Get matches
    matches = mol.GetSubstructMatches(hemiacetal_pattern)
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    
    # Patterns to exclude
    glycoside_pattern = Chem.MolFromSmarts("[OR0][C;R]1[O;R][C;R][C;R][C;R][C;R]1")
    lactone_pattern = Chem.MolFromSmarts("[O;R][C;R](=O)")
    
    valid_lactol = False
    for match in matches:
        # Check first two atoms (O and C) are in same ring
        rings_with_o = set(i for i, ring in enumerate(ring_info.AtomRings()) 
                          if match[0] in ring)
        rings_with_c = set(i for i, ring in enumerate(ring_info.AtomRings()) 
                          if match[1] in ring)
        
        common_rings = rings_with_o.intersection(rings_with_c)
        if not common_rings:
            continue
            
        # Get the smallest ring containing both atoms
        ring_size = min(len(ring_info.AtomRings()[i]) for i in common_rings)
        
        # Typical lactol rings are 5-7 membered
        if ring_size < 4 or ring_size > 8:
            continue
            
        # Check if carbon is sp3 hybridized
        c_atom = mol.GetAtomWithIdx(match[1])
        if c_atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
            
        valid_lactol = True
        break

    if not valid_lactol:
        return False, "No valid lactol ring found"

    # Exclude molecules that are primarily glycosides
    if mol.HasSubstructMatch(glycoside_pattern):
        # Count oxygen atoms
        o_count = len([a for a in mol.GetAtoms() if a.GetSymbol() == 'O'])
        # If molecule has many oxygens and glycoside pattern, it's likely a glycoside
        if o_count > 5:
            return False, "Structure appears to be a glycoside rather than a lactol"
    
    # Exclude molecules that are primarily lactones
    if mol.HasSubstructMatch(lactone_pattern):
        lactone_matches = len(mol.GetSubstructMatches(lactone_pattern))
        hemiacetal_matches = len(matches)
        if lactone_matches >= hemiacetal_matches:
            return False, "Structure appears to be a lactone rather than a lactol"

    # Additional check for complex polycyclic systems
    if len(ring_info.AtomRings()) > 5:
        # For complex systems, require clear hemiacetal character
        if len(matches) < 2:  # Require multiple hemiacetal groups
            return False, "Complex ring system without clear lactol character"

    return True, "Contains cyclic hemiacetal structure characteristic of lactols"