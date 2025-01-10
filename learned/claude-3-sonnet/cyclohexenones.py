"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: cyclohexenones
Definition: Any six-membered alicyclic ketone having one double bond in the ring.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import RingInfo

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule contains a cyclohexenone structure.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_cyclohexenone, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific SMARTS patterns for cyclohexenone
    patterns = [
        # Basic cyclohexenone pattern with explicit non-aromatic atoms
        "[C;!a]1[C;!a][C;!a]=[C;!a][C;!a][C;!a]1(=O)",
        # Alternative pattern with double bond in different position
        "[C;!a]1[C;!a]=[C;!a][C;!a][C;!a][C;!a]1(=O)",
        # Pattern for bridged systems
        "[C;!a]1[C;!a][C;!a]=[C;!a][C;!a]([C;!a]1=O)",
        # Pattern for substituted systems
        "[C;!a]1([*,H])[C;!a]([*,H])[C;!a]([*,H])=[C;!a]([*,H])[C;!a]([*,H])[C;!a]1(=O)",
    ]

    ring_info = mol.GetRingInfo()
    
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is None:
            continue
            
        matches = mol.GetSubstructMatches(patt)
        for match in matches:
            ring_atoms = set(match)
            
            # Skip if any atom is aromatic
            if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring_atoms):
                continue
                
            # Get all rings containing these atoms
            rings_containing_match = []
            for ring in ring_info.AtomRings():
                if any(idx in ring for idx in match):
                    rings_containing_match.append(set(ring))
            
            # Check if this is part of a larger conjugated system
            conjugated = False
            for ring in rings_containing_match:
                if len(ring.intersection(ring_atoms)) > 0 and len(ring) != 6:
                    conjugated = True
                    break
            if conjugated:
                continue

            # Count and validate bonds in the ring
            double_bonds = 0
            ketone_bonds = 0
            other_bonds = 0
            
            for bond in mol.GetBonds():
                start_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                
                # Only consider bonds where both atoms are in our ring
                if start_idx in ring_atoms and end_idx in ring_atoms:
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        start_atom = mol.GetAtomWithIdx(start_idx)
                        end_atom = mol.GetAtomWithIdx(end_idx)
                        
                        # Check if it's a ketone bond
                        if start_atom.GetAtomicNum() == 8 or end_atom.GetAtomicNum() == 8:
                            ketone_bonds += 1
                        else:
                            double_bonds += 1
                    else:
                        other_bonds += 1

            # Validate bond counts
            if (double_bonds == 1 and ketone_bonds == 1 and other_bonds == 4):
                # Additional check for proper ketone position relative to double bond
                ketone_atom = None
                double_bond_atoms = set()
                
                for atom_idx in ring_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if any(n.GetAtomicNum() == 8 for n in atom.GetNeighbors()):
                        ketone_atom = atom_idx
                    for bond in atom.GetBonds():
                        if (bond.GetBondType() == Chem.BondType.DOUBLE and 
                            bond.GetBeginAtomIdx() in ring_atoms and 
                            bond.GetEndAtomIdx() in ring_atoms and
                            bond.GetBeginAtom().GetAtomicNum() == 6 and 
                            bond.GetEndAtom().GetAtomicNum() == 6):
                            double_bond_atoms.add(bond.GetBeginAtomIdx())
                            double_bond_atoms.add(bond.GetEndAtomIdx())
                
                if ketone_atom is not None and double_bond_atoms:
                    return True, "Contains cyclohexenone ring (6-membered ring with one ketone and one double bond)"

    return False, "No valid cyclohexenone pattern found"