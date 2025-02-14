"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:36973 disaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdchem
from typing import Tuple

def is_disaccharide(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound where two monosaccharides are joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for glycosidic bond (-O-)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond (-O-) found"
    
    # Identify ring systems
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Count monosaccharide rings (5 or 6 members)
    monosaccharide_rings = [ring for ring in rings if len(ring) in [5, 6]]
    n_monosaccharides = len(monosaccharide_rings)
    
    if n_monosaccharides != 2:
        return False, f"Found {n_monosaccharides} monosaccharide rings, expected 2"
    
    # Check if the two monosaccharide rings are connected by a glycosidic bond
    for ring1, ring2 in ((ring1, ring2) for i, ring1 in enumerate(monosaccharide_rings)
                         for ring2 in monosaccharide_rings[i+1:]):
        shared_atoms = set(ring1) & set(ring2)
        if len(shared_atoms) == 1:
            shared_atom = mol.GetAtomWithIdx(list(shared_atoms)[0])
            if shared_atom.GetAtomicNum() == 8:  # Oxygen
                # Check for additional structural features if needed
                return True, "Contains two monosaccharide rings connected by a glycosidic bond"
    
    return False, "Monosaccharide rings not connected by a glycosidic bond"