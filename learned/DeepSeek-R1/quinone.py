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
    
    # Check for conjugated cyclic dione structure using multiple SMARTS patterns
    patterns = [
        # Benzoquinone-like (6-membered ring)
        Chem.MolFromSmarts('[O]=C1C=CC(=O)C=C1'),
        # Anthraquinone-like (fused rings)
        Chem.MolFromSmarts('[O]=C1C2=C(C(=O)C=C1)C=CC=C2'),
        # Naphthoquinone-like (two fused rings)
        Chem.MolFromSmarts('[O]=C1C(=O)C2=CC=CC=C2C=C1'),
        # General conjugated cyclic dione (any ring size)
        Chem.MolFromSmarts('[O]=C1~*~*~*~*~*1'),  # Adjust as needed
    ]
    
    for pattern in patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return True, "Contains conjugated cyclic dione structure"
    
    # Additional check for conjugated systems not captured by SMARTS
    # Get all rings containing at least two ketones
    rings = mol.GetRingInfo().AtomRings()
    ketone_atoms = {match.GetAtomIdx() for match in mol.GetSubstructMatches(Chem.MolFromSmarts('[#6]=O'))}
    for ring in rings:
        ring_ketones = [atom for atom in ring if atom in ketone_atoms]
        if len(ring_ketones) >= 2:
            # Check if the path between any two ketones in the ring is conjugated
            for i in range(len(ring_ketones)):
                for j in range(i+1, len(ring_ketones)):
                    path = Chem.GetShortestPath(mol, ring_ketones[i], ring_ketones[j])
                    conjugated = True
                    for k in range(len(path)-1):
                        bond = mol.GetBondBetweenAtoms(path[k], path[k+1])
                        if not bond.GetIsConjugated():
                            conjugated = False
                            break
                    if conjugated:
                        return True, "Contains conjugated cyclic dione structure"
    
    return False, "No conjugated cyclic dione structure found"