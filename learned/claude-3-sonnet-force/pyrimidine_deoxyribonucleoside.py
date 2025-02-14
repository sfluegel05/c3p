"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: CHEBI:35808 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside is a deoxyribonucleoside containing a pyrimidine base.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyrimidine rings
    pyrimidine_rings = [ring for ring in AllChem.GetSymmSSSR(mol) if sum(1 for atom in ring if atom.GetAtomicNum() in [7, 6]) == 5]
    if not pyrimidine_rings:
        return False, "No pyrimidine ring found"
    
    # Check for deoxyribose sugar
    def is_deoxyribose(ring):
        if not ring.GetNumAtoms() == 5:
            return False
        
        # Check for oxygen and carbon atoms
        oxygen_count = sum(1 for atom in ring if atom.GetAtomicNum() == 8)
        carbon_count = sum(1 for atom in ring if atom.GetAtomicNum() == 6)
        if oxygen_count != 1 or carbon_count != 4:
            return False
        
        # Check for attached hydroxyl groups and CH2OH group
        hydroxyl_count = sum(1 for atom in ring if atom.GetTotalNumHs() == 1 and atom.GetAtomicNum() == 8)
        ch2oh_count = sum(1 for atom in ring if atom.GetDegree() == 3 and atom.GetAtomicNum() == 6 and sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1) == 1)
        if not (hydroxyl_count == 2 and ch2oh_count == 1):
            return False
        
        return True
    
    deoxyribose_rings = [ring for ring in AllChem.GetSymmSSSR(mol) if is_deoxyribose(ring)]
    if not deoxyribose_rings:
        return False, "No deoxyribose sugar found"
    
    # Check if pyrimidine and deoxyribose are connected
    for pyrimidine_ring in pyrimidine_rings:
        for deoxyribose_ring in deoxyribose_rings:
            pyrimidine_atoms = set(pyrimidine_ring)
            deoxyribose_atoms = set(deoxyribose_ring)
            if any(atom_set & deoxyribose_atoms for atom_set in [set(atom.GetNeighbors()) for atom in pyrimidine_atoms]):
                return True, "Contains a pyrimidine base connected to a deoxyribose sugar"
    
    return False, "Pyrimidine and deoxyribose not connected"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35808',
        'name': 'pyrimidine deoxyribonucleoside',
        'definition': 'A deoxyribonucleoside containing a pyrimidine base.'
    },
    # ... (rest of the metadata)
}