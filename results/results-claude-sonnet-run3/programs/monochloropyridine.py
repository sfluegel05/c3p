from rdkit import Chem
from rdkit.Chem import AllChem

def is_monochloropyridine(smiles: str):
    """
    Determines if a molecule contains exactly one chlorine atom attached to a pyridine ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monochloropyridine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Find all pyridine rings
    pyridine_rings = []
    for atom in mol.GetAtoms():
        # Check if atom is nitrogen in aromatic ring
        if atom.GetIsAromatic() and atom.GetSymbol() == 'N':
            # Get the ring atoms
            ring = []
            for bond in atom.GetBonds():
                if bond.GetIsAromatic():
                    other_atom = bond.GetOtherAtom(atom)
                    ring.append(other_atom.GetIdx())
                    for next_bond in other_atom.GetBonds():
                        if next_bond.GetIsAromatic():
                            next_atom = next_bond.GetOtherAtom(other_atom)
                            if next_atom.GetIdx() not in ring:
                                ring.append(next_atom.GetIdx())
            if len(ring) == 5:  # 5 carbons + 1 nitrogen = 6-membered ring
                ring.append(atom.GetIdx())
                pyridine_rings.append(ring)

    if not pyridine_rings:
        return False, "No pyridine rings found"

    # Check each pyridine ring for chlorine substituents
    chloro_pyridines = []
    for ring in pyridine_rings:
        chlorine_count = 0
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        
        for ring_atom in ring_atoms:
            for neighbor in ring_atom.GetNeighbors():
                if neighbor.GetSymbol() == 'Cl' and neighbor.GetIdx() not in ring:
                    chlorine_count += 1
                    
        if chlorine_count == 1:
            chloro_pyridines.append(ring)
            
    if not chloro_pyridines:
        return False, "No pyridine rings with exactly one chlorine substituent found"
        
    return True, f"Found {len(chloro_pyridines)} pyridine ring(s) with exactly one chlorine substituent"
# Pr=1.0
# Recall=0.9411764705882353