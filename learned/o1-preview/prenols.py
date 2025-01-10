"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH
    in which the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Allow atoms commonly found in prenols
    allowed_atoms = {1, 6, 7, 8, 15}  # H, C, N, O, P
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains atom not typical of prenols: {atom.GetSymbol()}"
    
    # Check for a hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group (-OH) found"

    # Define isoprene unit SMARTS pattern
    isoprene_pattern = Chem.MolFromSmarts("C(=C(C))CC")
    
    # Find all isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)
    if num_isoprene_units == 0:
        return False, "No isoprene units found"
    
    # Check if isoprene units are connected head-to-tail
    # Get the atom indices of all isoprene units
    isoprene_atoms = set()
    for match in isoprene_matches:
        isoprene_atoms.update(match)
    
    # Generate a subgraph containing only isoprene atoms
    isoprene_submol = Chem.PathToSubmol(mol, list(isoprene_atoms))
    
    # Check that isoprene units form a continuous chain
    # Count the number of connected components in the isoprene subgraph
    num_components = Chem.GetMolFrags(isoprene_submol, asMols=False, sanitizeFrags=False)
    if len(num_components) > 1:
        return False, "Isoprene units are not connected head-to-tail"
    
    # Check for terminal -OH group attached to the isoprene chain
    # Find hydroxyl groups in the molecule
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    # Check if any hydroxyl group is attached to the isoprene chain
    isoprene_neighbors = set()
    for atom_idx in isoprene_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            isoprene_neighbors.add(neighbor.GetIdx())
    hydroxyl_attached = False
    for match in hydroxyl_matches:
        hydroxyl_atom_idx = match[0]
        if hydroxyl_atom_idx in isoprene_neighbors:
            hydroxyl_attached = True
            break
    if not hydroxyl_attached:
        return False, "No terminal hydroxyl group attached to isoprene chain"
    
    # All checks passed
    return True, f"Molecule contains {num_isoprene_units} isoprene units with terminal hydroxyl group"
    

__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'prenols',
        'definition': "Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).",
        'parents': []
    }
}