"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:32854 secondary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains a secondary amine group.
    A secondary amine has exactly one N-H bond and two non-hydrogen substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # SMARTS patterns to exclude
    exclude_patterns = [
        "[nX2]",  # aromatic nitrogen in 5-membered ring
        "[NX3]C(=O)",  # amide
        "[N+](=O)[O-]",  # nitro
        "[N]=[O]",  # nitroso
        "[NX2]=[C,N]",  # imine
        "[N-][N+]#N",  # azide
        "n1cccc1",  # nitrogen in pyrrole
        "[n;H0]",  # aromatic N without H
        "[NH0]",  # N without any H
        "[N+]",  # quaternary N
        "[N-]",  # negatively charged N
    ]

    # Create molecules from SMARTS patterns
    exclude_mols = [Chem.MolFromSmarts(pattern) for pattern in exclude_patterns]

    # Secondary amine pattern: N with exactly one H and at least 2 non-H neighbors
    sec_amine_pattern = Chem.MolFromSmarts("[NX3H1][!H]")
    
    if not mol.HasSubstructMatch(sec_amine_pattern):
        return False, "No nitrogen with one hydrogen and two substituents found"

    # Get matches for secondary amine pattern
    matches = mol.GetSubstructMatches(sec_amine_pattern)
    
    for match in matches:
        N_idx = match[0]
        N_atom = mol.GetAtomWithIdx(N_idx)
        
        # Skip if nitrogen is part of any excluded pattern
        if any(mol.HasSubstructMatch(pattern) for pattern in exclude_mols):
            continue
            
        # Count non-hydrogen neighbors
        heavy_neighbors = len([n for n in N_atom.GetNeighbors() 
                             if n.GetAtomicNum() != 1])
        
        # Secondary amine must have exactly 2 non-H substituents
        if heavy_neighbors == 2:
            # Check that the nitrogen is not part of a ring system
            if not N_atom.IsInRing() or N_atom.IsInRingSize(6):  # Allow 6-membered rings
                return True, "Contains secondary amine group (NH with two non-H substituents)"
            
    return False, "No valid secondary amine group found"