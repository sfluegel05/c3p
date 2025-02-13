"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18020 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondType

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with an aldehyde group at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose), existing as a linear or cyclic structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 6 carbon atoms and sufficient oxygen atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count != 6 or o_count < 5:
        return False, "Does not contain 6 carbon atoms or sufficient oxygen atoms for a hexose"
    
    # Define SMARTS patterns for aldohexose and ketohexose
    aldohexose_pattern = Chem.MolFromSmarts("[CH2](C(=O)C[CH][CH][CH][CH]O)")  # Linear aldohexose
    ketohexose_pattern = Chem.MolFromSmarts("[CH](C(=O)C[CH][CH][CH]O)")  # Linear ketohexose
    pyranose_pattern = Chem.MolFromSmarts("O1C[CH][CH][CH][CH]C1")  # Pyranose ring
    furanose_pattern = Chem.MolFromSmarts("O1C[CH][CH]C[CH]1")  # Furanose ring
    
    # Check for linear aldohexose or ketohexose
    if mol.HasSubstructMatch(aldohexose_pattern):
        return True, "Linear aldohexose structure found"
    if mol.HasSubstructMatch(ketohexose_pattern):
        return True, "Linear ketohexose structure found"
    
    # Check for cyclic pyranose or furanose structures
    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        # Check for aldehyde or ketone group in ring
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon atom
                bond_types = [bond.GetBondType() for bond in atom.GetBonds()]
                if BondType.DOUBLE in bond_types:
                    # Check if connected to oxygen (aldehyde or ketone)
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 8:
                            # Check position of aldehyde or ketone group
                            ring_atoms = mol.GetRingInfo().AtomRings[mol.GetRingInfo().IsCyclic(atom.GetIdx())][1]
                            if len(ring_atoms) == 6:
                                ring_pos = ring_atoms.index(atom.GetIdx())
                                if ring_pos == 0:
                                    return True, "Aldopyranose or aldofuranose structure found"
                                elif ring_pos == 1:
                                    return True, "Ketopyranose or ketofuranose structure found"
    
    return False, "No hexose structure found"