"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide
Definition: Any aromatic amide obtained by acylation of aniline.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide has an N-phenyl amide group (Ph-NH-C(=O)-R).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for N-phenyl amide:
    # [c] - aromatic carbon
    # [NH] - nitrogen with one hydrogen
    # [CX3](=[OX1]) - carbon with 3 connections, double bonded to oxygen
    anilide_pattern = Chem.MolFromSmarts('[c][NH][CX3](=[OX1])')
    
    # Find all matches
    matches = mol.GetSubstructMatches(anilide_pattern)
    
    if not matches:
        return False, "No N-phenyl amide group found"
    
    # For each match, verify it's part of a proper benzene ring
    valid_matches = 0
    for match in matches:
        aromatic_carbon = match[0]
        
        # Get the ring that contains this carbon
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()
        
        # Check if the aromatic carbon is part of a 6-membered ring
        for ring in rings:
            if aromatic_carbon in ring and len(ring) == 6:
                # Verify all atoms in the ring are aromatic carbons
                ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
                if all(atom.GetIsAromatic() and atom.GetAtomicNum() == 6 for atom in ring_atoms):
                    valid_matches += 1
                    break
    
    if valid_matches == 0:
        return False, "No N-phenyl amide group found with proper benzene ring"
    
    return True, f"Found {valid_matches} N-phenyl amide group(s) with proper benzene ring(s)"