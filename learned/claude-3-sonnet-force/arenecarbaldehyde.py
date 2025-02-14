"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: CHEBI:35815 arenecarbaldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find aromatic rings
    aromatic_rings = mol.GetRingInfo().AromaticRings()
    
    # Check if any aldehyde group is attached to an aromatic ring
    aldehyde_attached_to_aromatic = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetDegree() == 3 and sum(1 for a in atom.GetNeighbors() if a.GetAtomicNum() == 8) == 1:
            # Atom is an aldehyde carbon
            neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() != 8]
            for neighbor in neighbors:
                # Check if neighbor is part of an aromatic ring, or if there's a short aliphatic chain
                # connecting the aldehyde group to an aromatic ring
                for ring in aromatic_rings:
                    if neighbor.IsInRingOfSize(len(ring)):
                        aldehyde_attached_to_aromatic = True
                        break
                    else:
                        # Check if neighbor is part of a short aliphatic chain leading to an aromatic ring
                        chain_atoms = [neighbor]
                        for _ in range(3):  # Maximum chain length of 3
                            next_atoms = []
                            for chain_atom in chain_atoms:
                                for neighbor in chain_atom.GetNeighbors():
                                    if neighbor.IsInRing() and any(neighbor.IsInRingOfSize(len(r)) for r in aromatic_rings):
                                        aldehyde_attached_to_aromatic = True
                                        break
                                    elif neighbor.GetAtomicNum() in (6, 8):  # Carbon or oxygen
                                        next_atoms.append(neighbor)
                            chain_atoms = next_atoms
                            if aldehyde_attached_to_aromatic:
                                break

    if aldehyde_attached_to_aromatic:
        return True, "Molecule contains an aldehyde group attached to an aromatic moiety"
    else:
        return False, "No aldehyde group attached to an aromatic moiety found"