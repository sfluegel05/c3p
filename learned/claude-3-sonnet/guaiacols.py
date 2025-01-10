"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule contains a guaiacol moiety based on its SMILES string.
    A guaiacol is a phenol with a methoxy group in the ortho position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains guaiacol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible SMARTS patterns for guaiacol core:
    # Match aromatic carbon with OH and adjacent carbon with OMe
    # Try both orientations (OH-OMe and OMe-OH)
    patterns = [
        # Pattern 1: OH then OMe
        '[c;$(c1ccccc1)]([OH1])[$(c([OC]))]-,=',
        # Pattern 2: OMe then OH
        '[c;$(c1ccccc1)]([OCH3])[$(c([OH1]))]-,='
    ]
    
    for pattern in patterns:
        guaiacol_pattern = Chem.MolFromSmarts(pattern)
        if guaiacol_pattern is None:
            continue
            
        matches = mol.GetSubstructMatches(guaiacol_pattern)
        if matches:
            # Verify each match
            for match in matches:
                # Get the matched atoms
                core_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
                
                # Check if the carbons are part of an aromatic ring
                ring_info = mol.GetRingInfo()
                if all(ring_info.IsAtomInRingOfSize(atom.GetIdx(), 6) for atom in core_atoms 
                      if atom.GetAtomicNum() == 6):
                    return True, "Contains phenol with ortho methoxy group (guaiacol moiety)"
    
    # Additional check with a simpler pattern
    simple_pattern = Chem.MolFromSmarts('c(O)c(OC)')
    if mol.HasSubstructMatch(simple_pattern):
        return True, "Contains phenol with ortho methoxy group (guaiacol moiety)"
                
    return False, "No guaiacol pattern found"