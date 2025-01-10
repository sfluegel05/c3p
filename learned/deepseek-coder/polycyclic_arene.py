"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: CHEBI:33853 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene is a polycyclic aromatic hydrocarbon, which may have substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check if there are at least 2 rings
    if len(rings) < 2:
        return False, "Not enough rings to be polycyclic"

    # Check if the molecule has a polycyclic aromatic core
    # Create a copy of the molecule and remove all non-carbon atoms
    core_mol = Chem.RWMol(mol)
    for atom in core_mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            core_mol.RemoveAtom(atom.GetIdx())
    
    # Check if the core has at least 2 aromatic rings
    core_ring_info = core_mol.GetRingInfo()
    core_rings = core_ring_info.AtomRings()
    
    aromatic_rings = 0
    for ring in core_rings:
        if all(core_mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring):
            aromatic_rings += 1
    
    if aromatic_rings < 2:
        return False, "Core structure does not have enough aromatic rings"

    # Check if the core is polycyclic (at least 2 rings share atoms)
    shared_atoms = 0
    for i in range(len(core_rings)):
        for j in range(i + 1, len(core_rings)):
            shared_atoms = len(set(core_rings[i]).intersection(core_rings[j]))
            if shared_atoms > 0:
                break
        if shared_atoms > 0:
            break
    else:
        return False, "Core rings are not fused"

    return True, "Contains a polycyclic aromatic hydrocarbon core"