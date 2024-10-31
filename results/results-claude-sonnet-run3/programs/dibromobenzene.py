from rdkit import Chem
from rdkit.Chem import AllChem

def is_dibromobenzene(smiles: str):
    """
    Determines if a molecule is a dibromobenzene (benzene with exactly two bromo substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dibromobenzene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find aromatic rings
    rings = mol.GetRingInfo()
    
    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check each aromatic ring for exactly two bromine substituents
    for ring in aromatic_rings:
        ring_atoms = set(ring)
        br_count = 0
        
        # Check each atom in the ring for bromine neighbors
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atoms and neighbor.GetSymbol() == 'Br':
                    br_count += 1
        
        # If we find a ring with exactly 2 bromine substituents, return True
        if br_count == 2:
            return True, "Benzene ring with exactly two bromine substituents found"
            
    return False, "No benzene ring with exactly two bromine substituents found"
# Pr=0.6666666666666666
# Recall=1.0