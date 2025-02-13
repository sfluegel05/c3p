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

    # Basic cyclic hemiacetal patterns
    patterns = [
        # Ring oxygen connected to carbon with OH or OR
        "[O;R][C;R]([OH1,OR])",
        
        # Ring oxygen connected to carbon with OH, more specific
        "[O;R][C;R]([OH1])",
        
        # Hemiacetal in equilibrium with ketone/aldehyde form
        "[O;R][C;R](=O)",
        
        # Common natural product lactol pattern
        "[O;R][C;R]([OH1,OR])[C,H]",
        
        # Bridged lactol pattern
        "[O;R][C;R]([OH1,OR])[C;R]",
        
        # Fused ring lactol pattern
        "[O;R][C;R]([OH1,OR])[C;R0]"
    ]

    # Convert patterns to RDKit SMARTS objects
    smarts_patterns = [Chem.MolFromSmarts(p) for p in patterns]

    # Check for matches
    matches = []
    for pattern in smarts_patterns:
        if pattern is not None:  # Ensure valid SMARTS pattern
            matches.extend(mol.GetSubstructMatches(pattern))

    if not matches:
        return False, "No lactol structure found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    
    # Verify that matched structures are part of rings
    valid_lactol = False
    for match in matches:
        # Check first two atoms (O and C) are in same ring
        rings_with_o = set(i for i, ring in enumerate(ring_info.AtomRings()) 
                          if match[0] in ring)
        rings_with_c = set(i for i, ring in enumerate(ring_info.AtomRings()) 
                          if match[1] in ring)
        
        if rings_with_o.intersection(rings_with_c):
            valid_lactol = True
            break

    if not valid_lactol:
        return False, "Matched atoms not in same ring"

    # Additional checks for specific cases
    # Pattern to identify simple glycosides
    glycoside = Chem.MolFromSmarts("[OR0][C;R]1[O;R][C;R][C;R][C;R][C;R]1")
    
    # If molecule only matches glycoside pattern and has typical sugar structure,
    # it's probably a glycoside rather than a lactol
    if (mol.HasSubstructMatch(glycoside) and 
        len(matches) == 1 and 
        any(a.GetSymbol() == 'O' for a in mol.GetAtoms()) > 5):
        return False, "Structure appears to be a glycoside rather than a lactol"

    return True, "Contains cyclic hemiacetal structure characteristic of lactols"