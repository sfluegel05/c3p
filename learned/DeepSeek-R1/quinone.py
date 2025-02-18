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
    A quinone has a fully conjugated cyclic dione structure derived from aromatic compounds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all ketone groups (C=O)
    ketone_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[#6]=O'))
    if len(ketone_matches) < 2:
        return False, f"Found {len(ketone_matches)} ketone groups, need at least 2"
    
    # Check for conjugated cyclic dione structure using SMARTS patterns
    patterns = [
        # Benzoquinone (6-membered)
        Chem.MolFromSmarts('[O]=C1C=CC(=O)C=C1'),
        # Naphthoquinone (two fused rings)
        Chem.MolFromSmarts('[O]=C1C(=O)C2=CC=CC=C2C=C1'),
        # Anthraquinone (three fused rings)
        Chem.MolFromSmarts('[O]=C1C2=C(C(=O)C=C1)C=CC=C2'),
        # General conjugated cyclic dione (any even-membered ring)
        Chem.MolFromSmarts('[O]=C1C(=O)C=CC=C1'),  # 6-membered variant
        Chem.MolFromSmarts('[O]=C1C(=O)CC=C1'),     # 5-membered (if exists)
    ]
    
    for pattern in patterns:
        if pattern and mol.HasSubstructMatch(pattern):
            return True, "Contains conjugated cyclic dione structure"
    
    # Additional check for rings with at least two ketones in conjugation
    rings = mol.GetRingInfo().AtomRings()
    ketone_atoms = {match[0] for match in ketone_matches}  # Fix: extract first atom from each match tuple
    
    for ring in rings:
        ring_ketones = [atom for atom in ring if atom in ketone_atoms]
        if len(ring_ketones) >= 2:
            # Check all pairs of ketones in the ring for conjugated path
            for i in range(len(ring_ketones)):
                for j in range(i+1, len(ring_ketones)):
                    a1, a2 = ring_ketones[i], ring_ketones[j]
                    path = Chem.GetShortestPath(mol, a1, a2)
                    if not path:
                        continue
                    conjugated = all(mol.GetBondBetweenAtoms(path[k], path[k+1]).GetIsConjugated() 
                                    for k in range(len(path)-1))
                    if conjugated:
                        return True, "Conjugated dione in cyclic system"
    
    return False, "No conjugated cyclic dione structure found"