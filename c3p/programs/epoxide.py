"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:52092 epoxide
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, GetSSSR

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a cyclic ether with a 3-membered ring containing an oxygen atom bonded to two carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        # Generate all possible rings (including non-SSSR)
        mol.UpdatePropertyCache(strict=False)
        Chem.GetSSSR(mol)
    except:
        return False, "Error processing molecule structure"
    
    rings = mol.GetRingInfo().AtomRings()
    
    for ring in rings:
        if len(ring) == 3:
            # Check each atom in the 3-membered ring
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 8:  # Oxygen atom
                    # Get neighboring atoms in the ring
                    neighbors = [mol.GetAtomWithIdx(i) for i in ring if i != atom_idx]
                    
                    # Verify both neighbors are carbons
                    if all(n.GetAtomicNum() == 6 for n in neighbors):
                        # Check oxygen has exactly two single bonds
                        bonds = atom.GetBonds()
                        if len(bonds) == 2 and all(bond.GetBondType() == Chem.BondType.SINGLE for bond in bonds):
                            return True, "Contains a 3-membered ether ring (epoxide)"
    
    return False, "No 3-membered ether ring found"