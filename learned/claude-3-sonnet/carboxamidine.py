"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:51381 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group (RC(=NR)NR2).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains carboxamidine group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific carboxamidine pattern:
    # [CX3] - sp2 carbon
    # [NX2] - sp2 nitrogen (double bond)
    # [NX3] - sp3 nitrogen (single bond)
    # !$([C]=[O,S,P]) - not carbonyl, thiocarbonyl, or phosphorus
    # !$(C=N[N,O,S]) - not hydrazone, oxime, or thiohydrazone
    # !$(C([NX3])[NX3][NX3]) - not guanidine
    pattern = Chem.MolFromSmarts("""
        [CX3;!$([C]=[O,S,P]);!$(C=N[N,O,S]);!$(C([NX3])[NX3][NX3])]
        (
            [#6,#1;!$(C=N),!$(C=O)]
        )
        (
            [NX3;!$(N(C=O));!$(N(N=*));!$(N(C=N)[#7,#8,#16])]
        )
        =[NX2;!$(N=CN=*);!$(N=CC=O);!$(N=CN([#7,#8,#16]))]
    """)

    # Find matches
    matches = mol.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No carboxamidine group found"
    
    # Additional checks for problematic cases
    for match in matches:
        carbon_idx = match[0]
        
        # Get all atoms connected to the carbon
        neighbors = mol.GetAtomWithIdx(carbon_idx).GetNeighbors()
        
        # Check if part of problematic cyclic structures
        ring_info = mol.GetRingInfo()
        if ring_info.IsAtomInRingOfSize(carbon_idx, 5):
            ring_atoms = set()
            for ring in ring_info.AtomRings():
                if carbon_idx in ring:
                    ring_atoms.update(ring)
            
            # Check if both nitrogens are in same 5-membered ring
            n_count = 0
            for n_idx in [match[2], match[3]]:  # indices of N atoms
                if n_idx in ring_atoms:
                    n_count += 1
            if n_count == 2:
                continue  # Skip this match as it's likely a cyclic imine
    
        # Valid carboxamidine found
        return True, f"Contains carboxamidine group"
    
    return False, "No valid carboxamidine group found"