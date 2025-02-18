"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a five-membered cyclic ester (lactone ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all ester groups in rings
    # SMARTS pattern for ester group in a 5-membered ring:
    # O connected to C=O in a 5-atom ring where both oxygen atoms are in the ring
    gamma_lactone_pattern = Chem.MolFromSmarts("[O;R][C;R](=O)[C;R][C;R][C;R]")
    
    # Check for matches with the 5-membered ring ester pattern
    matches = mol.GetSubstructMatches(gamma_lactone_pattern)
    
    if not matches:
        return False, "No five-membered lactone ring found"
    
    # Verify ring size for each matched oxygen atom
    for match in matches:
        ring_info = mol.GetRingInfo()
        for atom_idx in match:
            # Get rings containing the matched atom
            rings = [r for r in ring_info.AtomRings() if atom_idx in r]
            for ring in rings:
                # Check if it's a 5-membered ring containing both ester oxygens
                if len(ring) == 5:
                    # Get atoms in the ring
                    ring_atoms = set(ring)
                    # Check if both ester oxygen and carbonyl oxygen are in this ring
                    ester_o = match[0]
                    carbonyl_o = [a.GetIdx() for a in mol.GetAtomWithIdx(match[1]).GetNeighbors() if a.GetAtomicNum() == 8][0]
                    if {ester_o, carbonyl_o}.issubset(ring_atoms):
                        return True, "Contains a five-membered lactone ring"
    
    return False, "No valid five-membered lactone ring found"