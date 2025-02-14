"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: CHEBI:36343 Cyclic fatty acid
A cyclic fatty acid is defined as any fatty acid containing anywhere in its structure a ring of atoms.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_atoms = mol.GetSubstructMatches(acid_pattern)
    if not acid_atoms:
        return False, "No carboxylic acid group found"

    # Identify rings (aromatic and aliphatic)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Check if the carboxylic acid group is part of a ring
    acid_in_ring = False
    for ring in rings:
        for atom_idx in acid_atoms[0]:
            if atom_idx in ring:
                acid_in_ring = True
                break

    # Check for carbon chains
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, "Too few carbons for a fatty acid"

    if acid_in_ring:
        return True, "Contains a carboxylic acid group as part of a ring structure"
    else:
        # Check if any ring is connected to the carboxylic acid group
        acid_connected_rings = []
        for ring in rings:
            for atom_idx in ring:
                if atom_idx in [x[0] for x in acid_atoms]:
                    acid_connected_rings.append(ring)
                    break

        if acid_connected_rings:
            return True, "Contains a carboxylic acid group and at least one ring (aromatic or aliphatic) connected to the fatty acid chain"
        else:
            return False, "No rings connected to or including the carboxylic acid group"

    # Additional checks or filters can be added here if needed