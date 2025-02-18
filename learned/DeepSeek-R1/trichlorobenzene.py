"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: CHEBI:36027 trichlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import BondType

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule is a trichlorobenzene based on its SMILES string.
    A trichlorobenzene has at least one benzene ring with exactly three chlorine substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a benzene ring with three chlorines, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ri = mol.GetRingInfo()
    benzene_rings = []
    
    # Iterate through all rings to find benzene rings (6-membered aromatic carbon rings)
    for ring in ri.AtomRings():
        if len(ring) != 6:
            continue  # Not a six-membered ring
        is_benzene = True
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check if atom is an aromatic carbon
            if not (atom.GetIsAromatic() and atom.GetSymbol() == 'C'):
                is_benzene = False
                break
        if is_benzene:
            benzene_rings.append(ring)
    
    if not benzene_rings:
        return False, "No benzene ring found"
    
    # Check each benzene ring for exactly three chlorine substituents
    for ring in benzene_rings:
        cl_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                # Ensure the bond is single (substituent, not part of a fused ring)
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() != BondType.SINGLE:
                    continue
                if neighbor.GetSymbol() == 'Cl':
                    cl_count += 1
        if cl_count == 3:
            return True, "Benzene ring with three chlorine substituents found"
    
    return False, "No benzene ring with three chlorines"