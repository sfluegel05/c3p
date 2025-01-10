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

    # Basic requirements - must have both oxygen and carbon
    if not all(mol.GetAtomWithIdx(i).GetAtomicNum() in [6, 8] for i in range(mol.GetNumAtoms())):
        return False, "Missing required elements (C, O)"

    # SMARTS pattern for guaiacol:
    # c1(O)c(OC)cccc1 - phenol with ortho methoxy
    # We'll try both orientations since they're equivalent
    pattern1 = Chem.MolFromSmarts('c1(O)c(OC)cccc1')
    pattern2 = Chem.MolFromSmarts('c1(OC)c(O)cccc1')
    
    if pattern1 is None or pattern2 is None:
        return None, "Error in SMARTS patterns"

    # Check for matches
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    
    all_matches = list(matches1) + list(matches2)
    
    if not all_matches:
        return False, "No guaiacol pattern found"

    # Verify each match
    for match in all_matches:
        # Get the atoms involved
        atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        
        # Verify aromaticity
        if not all(atom.GetIsAromatic() for atom in atoms if atom.GetAtomicNum() == 6):
            continue
            
        # Find the oxygen atoms
        o_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 8]
        
        # Must have exactly two oxygens
        if len(o_atoms) != 2:
            continue
            
        # Verify one is a hydroxyl (OH) and one is part of methoxy (OCH3)
        has_hydroxyl = False
        has_methoxy = False
        
        for o_atom in o_atoms:
            h_count = o_atom.GetTotalNumHs()
            c_neighbors = [n for n in o_atom.GetNeighbors() if n.GetAtomicNum() == 6]
            
            if h_count == 1 and len(c_neighbors) == 1:
                has_hydroxyl = True
            elif h_count == 0 and len(c_neighbors) == 2:
                # One neighbor should be aromatic, one should be sp3
                if (any(n.GetIsAromatic() for n in c_neighbors) and 
                    any(n.GetHybridization() == Chem.HybridizationType.SP3 for n in c_neighbors)):
                    has_methoxy = True
        
        if has_hydroxyl and has_methoxy:
            return True, "Contains phenol with ortho methoxy group (guaiacol moiety)"
            
    return False, "Contains required elements but not in correct configuration"