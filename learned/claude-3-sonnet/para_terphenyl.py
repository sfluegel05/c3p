"""
Classifies: CHEBI:75874 para-terphenyl
"""
"""
Classifies: para-terphenyl compounds
A ring assembly based on a 1,4-diphenylbenzene skeleton and its substituted derivatives.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_para_terphenyl(smiles: str):
    """
    Determines if a molecule is a para-terphenyl based on its SMILES string.
    Para-terphenyls have a central benzene ring with two phenyl groups attached in para positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a para-terphenyl, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for presence of exactly 3 benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    if len(benzene_matches) < 3:
        return False, "Less than 3 benzene rings found"
    
    # Look for para-terphenyl core structure
    # This pattern enforces:
    # - Central benzene ring
    # - Two phenyl groups in para positions
    # - Allows substitutions on any position
    para_terphenyl_pattern = Chem.MolFromSmarts("[cR1]1[cR1][cR1](-[cR1]2[cR1][cR1][cR1][cR1][cR1]2)[cR1][cR1](-[cR1]3[cR1][cR1][cR1][cR1][cR1]3)[cR1]1")
    
    if not mol.HasSubstructMatch(para_terphenyl_pattern):
        return False, "No para-terphenyl core structure found"
    
    # Get all rings
    ring_info = mol.GetRingInfo()
    rings = list(ring_info.AtomRings())
    
    # Check that we have exactly 3 six-membered aromatic rings
    aromatic_rings = [ring for ring in rings 
                     if len(ring) == 6 and 
                     all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
    
    if len(aromatic_rings) != 3:
        return False, f"Found {len(aromatic_rings)} aromatic 6-membered rings, need exactly 3"

    # Check that rings are not fused (should share no atoms)
    for i, ring1 in enumerate(aromatic_rings):
        for ring2 in aromatic_rings[i+1:]:
            shared_atoms = set(ring1).intersection(set(ring2))
            if len(shared_atoms) > 0:
                return False, "Contains fused rings"
    
    # Find the central ring (should be connected to both other rings)
    for ring in aromatic_rings:
        connections = 0
        ring_atoms = set(ring)
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms and neighbor.GetIsAromatic():
                    connections += 1
        if connections == 2:  # This is the central ring
            break
    else:
        return False, "Could not identify central ring with correct connectivity"

    return True, "Contains para-terphenyl core structure with appropriate substitution pattern"