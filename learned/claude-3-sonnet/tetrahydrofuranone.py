"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: tetrahydrofuranone - any oxolane having an oxo- substituent at any position 
on the tetrahydrofuran ring
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone contains a tetrahydrofuran ring with a ketone group directly
    attached to the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core patterns for tetrahydrofuranones
    patterns = [
        # Basic tetrahydrofuranone with ketone in ring
        "[O;R1]1[C;R1]([C;R1,c;R1])[C;R1,c;R1][C;R1,c;R1][C;R1]1(=O)",
        
        # Anhydride-type patterns (two ketones)
        "[O;R1]1[C;R1](=O)[C;R1,c;R1][C;R1,c;R1][C;R1]1(=O)",
        
        # Lactone patterns
        "[O;R1]1[C;R1](=O)[C;R1,c;R1][C;R1,c;R1][C;R1,c;R1]1",
        "[O;R1]1[C;R1,c;R1][C;R1](=O)[C;R1,c;R1][C;R1,c;R1]1",
        "[O;R1]1[C;R1,c;R1][C;R1,c;R1][C;R1](=O)[C;R1,c;R1]1",
        "[O;R1]1[C;R1,c;R1][C;R1,c;R1][C;R1,c;R1][C;R1]1(=O)"
    ]

    # Check each pattern
    for pattern in patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is None:
            continue
        
        matches = mol.GetSubstructMatches(pattern_mol)
        if matches:
            # Verify that we have a proper 5-membered ring with oxygen
            for match in matches:
                ring_atoms = set(match)
                ring_info = mol.GetRingInfo()
                
                # Check if these atoms form a 5-membered ring
                for ring in ring_info.AtomRings():
                    ring = set(ring)
                    if len(ring & ring_atoms) == 5:
                        # Verify presence of ketone group
                        for atom_idx in ring:
                            atom = mol.GetAtomWithIdx(atom_idx)
                            for bond in atom.GetBonds():
                                if bond.GetBondType() == Chem.BondType.DOUBLE:
                                    other_atom = bond.GetOtherAtom(atom)
                                    if other_atom.GetAtomicNum() == 8:  # Oxygen
                                        return True, "Contains tetrahydrofuran ring with ketone group"

    return False, "No tetrahydrofuranone core structure found"