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
    # - Must have benzene ring (c1ccccc1)
    # - Must have acetate group (-OC(=O)CH3) directly attached to the benzene ring
    # - The oxygen must be connected to exactly one aromatic carbon and one acetyl group
    phenyl_acetate_pattern = Chem.MolFromSmarts("c1c([OX2][CX3](=[OX1])[CH3])cccc1")
    
    # Negative patterns - structures that should not be counted
    negative_patterns = [
        Chem.MolFromSmarts("O=C1Oc2ccccc2C1"),  # coumarin core
        Chem.MolFromSmarts("C1=NCc2ccccc12"),    # indole core
        Chem.MolFromSmarts("O=C1C=Cc2ccccc2O1"), # benzofuranone
        Chem.MolFromSmarts("O=C1C=COc2ccccc21"), # benzopyranone
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

        # Verify the aromatic carbon is part of a benzene ring
        ring_info = mol.GetRingInfo()
        ring_size = 0
        for ring in ring_info.AtomRings():
            if aromatic_carbon.GetIdx() in ring:
                ring_size = len(ring)
                break
                
        if ring_size != 6:
            continue

        # Count number of aromatic atoms in the ring to ensure it's benzene
        ring_aromatic_count = sum(1 for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetIsAromatic())
        if ring_aromatic_count != 6:
            continue

        valid_matches.append(match)

    if not valid_matches:
        return False, "No valid phenyl acetate group found"

    # Success case
    if len(valid_matches) >= 1:
        return True, "Contains phenyl ring with acetate group properly attached to phenolic oxygen"

    return False, "Structure does not match phenyl acetate definition"