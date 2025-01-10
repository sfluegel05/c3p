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
    
    # SMARTS patterns for different types of beta-lactams
    patterns = [
        # Basic beta-lactam core (more general)
        ("[N]1[C][C](=[O])[C]1", "basic beta-lactam"),
        # Penicillin-type structure
        ("[N]1[C]2[C](=[O])[C]1[S][C]([C])([C])[CH]2", "penicillin-type"),
        # Cephalosporin-type structure
        ("[N]1[C]2[C](=[O])[C]1[S][C][C]2", "cephalosporin-type"),
        # Carbapenem-type structure
        ("[N]1[C]2[C](=[O])[C]1[C][C]2", "carbapenem-type"),
        # Monobactam-type structure
        ("[N]1[C]([C])[C](=[O])[C]1[S](=[O])(=[O])[O]", "monobactam-type")
    ]
    
    matches = []
    match_types = set()
    
    for smarts, desc in patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matches.extend(mol.GetSubstructMatches(pattern))
            match_types.add(desc)
    
    if not matches:
        return False, "No beta-lactam ring found"
    
    # Verify each match contains an amide group in a 4-membered ring
    valid_rings = 0
    ring_info = mol.GetRingInfo()
    
    for match in matches:
        # Get all rings containing these atoms
        ring_atoms = set(match)
        for ring in ring_info.AtomRings():
            ring_set = set(ring)
            if len(ring_set.intersection(ring_atoms)) >= 4:
                # Check for amide bond
                for atom_idx in ring:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomicNum() == 7:  # Nitrogen
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 6:  # Carbon
                                for bond in neighbor.GetBonds():
                                    other_atom = bond.GetOtherAtom(neighbor)
                                    if (other_atom.GetAtomicNum() == 8 and 
                                        bond.GetBondType() == Chem.BondType.DOUBLE):
                                        valid_rings += 1
                                        break
    
    if valid_rings == 0:
        return False, "Found potential ring but missing required amide bond"
    
    # Additional structural validation
    if len(match_types) > 0:
        match_desc = ", ".join(match_types)
        return True, f"Contains beta-lactam structure(s): {match_desc}"
    
    return True, "Contains beta-lactam ring with amide bond"