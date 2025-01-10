"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:32952 primary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule contains a primary amine group based on its SMILES string.
    A primary amine has an NH2 group connected to a carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and handle salts by splitting on dots
    mol_parts = smiles.split('.')
    for part in mol_parts:
        mol = Chem.MolFromSmiles(part)
        if mol is None:
            continue
            
        # Add explicit hydrogens
        mol = Chem.AddHs(mol)

        # Primary amine patterns:
        # [NX3H2][#6] matches NH2 connected to any carbon
        # Exclude N that's part of amide, but check carefully
        primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][#6]")
        amide_pattern = Chem.MolFromSmarts("[NX3H2][CX3](=[OX1])")
        
        if not mol.HasSubstructMatch(primary_amine_pattern):
            continue
            
        # Get matches
        amine_matches = mol.GetSubstructMatches(primary_amine_pattern)
        amide_matches = mol.GetSubstructMatches(amide_pattern)
        
        # Convert matches to sets of nitrogen atoms
        amine_nitrogens = {match[0] for match in amine_matches}
        amide_nitrogens = {match[0] for match in amide_matches}
        
        # Remove amide nitrogens
        true_primary_amines = amine_nitrogens - amide_nitrogens
        
        if not true_primary_amines:
            continue
            
        # Verify each potential primary amine
        for n_idx in true_primary_amines:
            n_atom = mol.GetAtomWithIdx(n_idx)
            
            # Must have exactly 2 hydrogens
            if n_atom.GetTotalNumHs() != 2:
                continue
                
            # Must have exactly one non-H neighbor
            heavy_neighbors = [neighbor for neighbor in n_atom.GetNeighbors() 
                             if neighbor.GetAtomicNum() != 1]
            if len(heavy_neighbors) != 1:
                continue
                
            # The neighbor must be carbon
            if heavy_neighbors[0].GetAtomicNum() != 6:
                continue
                
            # Check that nitrogen is not part of a special group like azide
            if n_atom.IsInRing():  # Optional: exclude ring nitrogens if needed
                ring_size = min(len(ring) for ring in mol.GetRingInfo().AtomRings() 
                              if n_idx in ring)
                if ring_size < 4:  # Exclude very small rings that might be azides
                    continue
                    
            # If we get here, we've found a valid primary amine
            return True, "Contains primary amine group(s)"
            
    return False, "No valid primary amine groups found"