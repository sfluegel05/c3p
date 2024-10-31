from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroxybenzaldehyde(smiles: str):
    """
    Determines if a molecule is a dihydroxybenzaldehyde (benzaldehyde with two hydroxy substituents).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroxybenzaldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Find benzene rings
    rings = mol.GetRingInfo()
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # For each aromatic ring, check if it has both an aldehyde and two hydroxy groups
    for ring_atoms in aromatic_rings:
        ring_atom_set = set(ring_atoms)
        hydroxy_count = 0
        has_aldehyde = False

        # Check each atom in the ring for substituents
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            
            # Check neighbors of ring atoms
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in ring_atom_set:
                    # Check for aldehyde group (C=O)
                    if neighbor.GetSymbol() == 'C' and neighbor.GetFormalCharge() == 0:
                        for n2 in neighbor.GetNeighbors():
                            if n2.GetSymbol() == 'O' and n2.GetDegree() == 1:
                                has_aldehyde = True
                                break
                    
                    # Check for hydroxy group (-OH)
                    if neighbor.GetSymbol() == 'O' and neighbor.GetDegree() == 1:
                        hydroxy_count += 1

        if has_aldehyde and hydroxy_count >= 2:
            return True, f"Found benzaldehyde with {hydroxy_count} hydroxy groups"

    return False, "No benzaldehyde with two hydroxy groups found"
# Pr=0.75
# Recall=1.0