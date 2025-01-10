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

    # First check for presence of at least 3 benzene rings
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    if len(benzene_matches) < 3:
        return False, "Less than 3 benzene rings found"
    
    # Check for para-terphenyl core structure
    # Central benzene with two phenyl groups in para positions
    # Note: The central ring can have substituents, so we use wildcards [*] for non-hydrogen atoms
    para_terphenyl_pattern = Chem.MolFromSmarts("[$(c1c([*])c([*])c(c([*])c1[*])c2ccccc2)]-c3ccccc3")
    
    if not mol.HasSubstructMatch(para_terphenyl_pattern):
        return False, "No para-terphenyl core structure found"
    
    # Additional check to ensure the phenyl groups are actually in para positions
    # by looking at the distance between the attachment points
    matches = mol.GetSubstructMatches(para_terphenyl_pattern)
    for match in matches:
        # Get the atoms of the central ring
        central_ring_atoms = set(match[:6])  # First 6 atoms should be the central ring
        
        # Find the carbons where the outer phenyl groups are attached
        attachment_points = []
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # If this atom is in the central ring and connected to a phenyl group
            if atom_idx in central_ring_atoms:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in central_ring_atoms and neighbor.GetIsAromatic():
                        attachment_points.append(atom_idx)
        
        # If we found exactly 2 attachment points and they are para to each other
        if len(attachment_points) == 2:
            # In a benzene ring, para positions are separated by 3 bonds
            path_length = len(Chem.GetShortestPath(mol, attachment_points[0], attachment_points[1])) - 1
            if path_length == 3:
                return True, "Contains para-terphenyl core structure with appropriate substitution pattern"
    
    return False, "Phenyl groups not in para positions"