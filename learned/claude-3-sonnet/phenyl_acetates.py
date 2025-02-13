"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: phenyl acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
           with the hydroxy group of any phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a phenyl acetate and the reason
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for phenyl acetate:
    # More specific pattern requiring phenolic oxygen
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1([OH0]-[CX3](=[OX1])[CH3])ccccc1")
    
    # Expanded negative patterns
    negative_patterns = [
        # Fused ring systems
        Chem.MolFromSmarts("O=C1Oc2ccccc2C1"),  # coumarin
        Chem.MolFromSmarts("C1=NCc2ccccc12"),    # indole
        Chem.MolFromSmarts("O=C1C=Cc2ccccc2O1"), # benzofuranone
        Chem.MolFromSmarts("O=C1C=COc2ccccc21"), # benzopyranone
        Chem.MolFromSmarts("c12ccccc1cccc2"),    # naphthalene
        Chem.MolFromSmarts("c12ccccc1occc2"),    # benzofuran
        Chem.MolFromSmarts("c12ccccc1Occc2"),    # benzodioxin
    ]

    # Check for matches
    matches = mol.GetSubstructMatches(phenyl_acetate_pattern)
    if not matches:
        return False, "No phenyl acetate group found"

    # Check for negative patterns
    for pattern in negative_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains excluded ring system"

    # For each match, verify it's a proper phenyl acetate
    valid_matches = []
    for match in matches:
        # Get the oxygen atom
        oxygen_idx = match[1]
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        
        # Get the aromatic carbon it's attached to
        aromatic_carbon = None
        for neighbor in oxygen.GetNeighbors():
            if neighbor.GetIsAromatic():
                aromatic_carbon = neighbor
                break
                
        if aromatic_carbon is None:
            continue

        # Verify the aromatic carbon is part of an isolated benzene ring
        ring_info = mol.GetRingInfo()
        ring_atoms = None
        for ring in ring_info.AtomRings():
            if aromatic_carbon.GetIdx() in ring:
                ring_atoms = ring
                break
                
        if ring_atoms is None or len(ring_atoms) != 6:
            continue

        # Verify ring is isolated (not part of fused system)
        is_isolated = True
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if any atom in ring is part of another ring
            other_rings = [r for r in ring_info.AtomRings() if atom_idx in r and r != ring_atoms]
            if other_rings:
                is_isolated = False
                break
                
        if not is_isolated:
            continue

        # Verify phenolic nature - oxygen should be attached to only aromatic carbon and acetyl
        if len(oxygen.GetNeighbors()) != 2:
            continue
            
        acetyl_carbon = None
        for neighbor in oxygen.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and not neighbor.GetIsAromatic():
                acetyl_carbon = neighbor
                break
                
        if acetyl_carbon is None or len(acetyl_carbon.GetNeighbors()) != 3:
            continue

        valid_matches.append(match)

    if not valid_matches:
        return False, "No valid phenyl acetate group found"

    return True, "Contains phenyl ring with acetate group properly attached to phenolic oxygen"