from rdkit import Chem
from rdkit.Chem import AllChem

def is_difluorobenzene(smiles: str):
    """
    Determines if a molecule contains a difluorobenzene moiety (benzene ring with exactly 2 fluorine substituents).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains difluorobenzene, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Find all benzene rings
    rings = mol.GetRingInfo()
    benzene_rings = []
    
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms) and \
               all(atom.GetSymbol() == 'C' for atom in atoms):
                benzene_rings.append(ring)
                
    if not benzene_rings:
        return False, "No benzene rings found"
        
    # For each benzene ring, check if it has exactly 2 fluorine substituents
    for ring in benzene_rings:
        ring_atoms = set(ring)
        f_count = 0
        
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms:
                    if neighbor.GetSymbol() == 'F':
                        f_count += 1
                        
        if f_count == 2:
            # Get the positions of the fluorines relative to each other
            f_positions = []
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in ring_atoms and neighbor.GetSymbol() == 'F':
                        f_positions.append(atom_idx)
            
            # Determine substitution pattern
            ring_list = list(ring)
            pos1 = ring_list.index(f_positions[0])
            pos2 = ring_list.index(f_positions[1])
            diff = abs(pos1 - pos2)
            
            if diff == 1 or (diff == 5 and (pos1 == 0 or pos2 == 0)):
                pattern = "1,2-difluoro"
            elif diff == 2 or diff == 4:
                pattern = "1,3-difluoro"
            else:
                pattern = "1,4-difluoro"
                
            return True, f"Contains {pattern}benzene moiety"
            
    return False, "No benzene ring with exactly 2 fluorine substituents found"
# Pr=1.0
# Recall=0.972972972972973