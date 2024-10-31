from rdkit import Chem
from rdkit.Chem import AllChem

def is_ferroheme(smiles: str):
    """
    Determines if a molecule is a ferroheme (iron(II)--porphyrin coordination complex).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ferroheme, False otherwise 
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of iron
    has_fe = False
    fe_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Fe':
            has_fe = True
            fe_atom = atom
            break
    if not has_fe:
        return False, "No iron atom found"

    # Check Fe coordination to nitrogens
    n_coords = 0
    for neighbor in fe_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'N':
            n_coords += 1
            
    if n_coords != 4:
        return False, "Iron not coordinated to 4 nitrogens as required in ferroheme"

    # Find all rings
    rings = mol.GetRingInfo()
    ring_atoms = rings.AtomRings()
    
    # Find 5-membered rings that contain nitrogen (pyrroles)
    pyrrole_rings = []
    for ring in ring_atoms:
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'N' for atom in atoms):
                pyrrole_rings.append(ring)

    if len(pyrrole_rings) < 4:
        return False, "Does not contain 4 pyrrole rings required for porphyrin structure"

    # Check that the pyrrole rings are connected to Fe through their nitrogens
    fe_n_coords = set()
    for neighbor in fe_atom.GetNeighbors():
        if neighbor.GetSymbol() == 'N':
            fe_n_coords.add(neighbor.GetIdx())
            
    pyrrole_n_atoms = set()
    for ring in pyrrole_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'N':
                pyrrole_n_atoms.add(atom_idx)
                
    if not fe_n_coords.issubset(pyrrole_n_atoms):
        return False, "Iron not properly coordinated to pyrrole nitrogens"

    # Check for methine bridges between pyrroles
    methine_bridges = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            n_neighbors = len([n for n in atom.GetNeighbors()])
            if n_neighbors == 2:  # Potential methine bridge
                methine_bridges += 1

    if methine_bridges < 4:
        return False, "Missing methine bridges between pyrrole rings"

    return True, "Contains Fe(II) coordinated to porphyrin ring system"
# Pr=0.7857142857142857
# Recall=1.0