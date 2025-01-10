"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: CHEBI:15903 D-hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines if a molecule is a D-hexose based on its SMILES string.
    A D-hexose is a hexose that has D-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check basic hexose requirements
    # Should have 6 carbons and at least 5 oxygens (including ring oxygen)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count != 6:
        return False, "Not exactly 6 carbons"
    if o_count < 5:
        return False, "Not enough oxygens for a hexose"

    # Generate 3D coordinates to better determine stereochemistry
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)

    # Find the carbon at position 5
    # In linear form: C5 is the 5th carbon from the aldehyde end
    # In cyclic form: C5 is the carbon before the ring oxygen
    position_5_carbon = None
    
    # Try to find linear form first
    aldehyde = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 1:  # Aldehyde carbon
            aldehyde = atom
            break
    
    if aldehyde:
        # Linear form - trace the chain
        current = aldehyde
        for i in range(5):
            neighbors = [n for n in current.GetNeighbors() if n.GetAtomicNum() == 6]
            if not neighbors:
                break
            current = neighbors[0]
            if i == 4:  # Position 5
                position_5_carbon = current
    else:
        # Cyclic form - find the carbon before the ring oxygen
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) == 5 or len(ring) == 6:  # Furanose or pyranose
                for i, atom_idx in enumerate(ring):
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetAtomicNum() == 8:  # Ring oxygen
                        # The previous atom in the ring is position 5
                        prev_idx = ring[i-1] if i > 0 else ring[-1]
                        position_5_carbon = mol.GetAtomWithIdx(prev_idx)
                        break

    if not position_5_carbon:
        return False, "Could not identify position 5 carbon"

    # Check stereochemistry at position 5
    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        if position_5_carbon.GetProp('_CIPCode') == 'R':
            return True, "Hexose with D-configuration at position 5"
    except:
        pass

    # Alternative check using chiral tag
    if position_5_carbon.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
        return True, "Hexose with D-configuration at position 5"

    return False, "No D-configuration found at position 5"