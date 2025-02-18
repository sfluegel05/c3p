"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    Quinones have a fully conjugated cyclic dione structure derived from aromatic compounds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all ketone groups (C=O not in ester/amide)
    ketone_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[#6]=[#8]'))
    if len(ketone_matches) != 2:
        return False, f"Found {len(ketone_matches)} ketone groups, need exactly 2"
    
    # Get rings containing the ketone atoms
    ketone_atoms = {match[0] for match in ketone_matches}
    rings = mol.GetRingInfo().AtomRings()
    candidate_rings = [ring for ring in rings if ketone_atoms.issubset(ring)]
    
    if not candidate_rings:
        return False, "Ketones not in the same ring"
    
    # Check conjugation between the two ketones in the ring
    k1, k2 = ketone_atoms
    path = Chem.GetShortestPath(mol, k1, k2)
    if not path:
        return False, "Ketones not connected"
    
    # Ensure all bonds in the path are conjugated
    conjugated = all(mol.GetBondBetweenAtoms(u, v).GetIsConjugated() 
                    for u, v in zip(path[:-1], path[1:]))
    if not conjugated:
        return False, "Ketones not conjugated"
    
    # Check if the ring is conjugated (aromatic or alternating bonds)
    ring = next(iter(candidate_rings))  # check first candidate ring
    ring_bonds = [mol.GetBondBetweenAtoms(u, v) for u, v in zip(ring, ring[1:] + ring[:1])]
    if not all(bond.GetIsConjugated() for bond in ring_bonds):
        return False, "Ring not fully conjugated"
    
    # Verify the two ketones are part of the conjugated ring system
    return True, "Contains conjugated cyclic dione structure"