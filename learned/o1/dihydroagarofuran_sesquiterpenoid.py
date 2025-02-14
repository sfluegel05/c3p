"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: dihydroagarofuran sesquiterpenoid
"""

from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    A dihydroagarofuran sesquiterpenoid has a fused tricyclic system consisting of two cyclopentane rings
    and a tetrahydrofuran ring (a five-membered ring containing an oxygen atom), forming the dihydroagarofuran skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sesquiterpenoid (15-carbon terpenoid)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 15:
        return False, f"Number of carbons ({num_carbons}) less than 15"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()

    # Initialize list to collect potential tetrahydrofuran rings (5-membered rings with oxygen)
    thf_rings = []  # List of tuples (ring index, atom indices, bond indices)

    # Collect all 5-membered rings containing an oxygen atom
    for ring_idx, ring_atoms in enumerate(atom_rings):
        if len(ring_atoms) == 5:
            # Check if ring contains an oxygen atom
            if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring_atoms):
                # Store the ring info
                thf_rings.append({
                    'idx': ring_idx,
                    'atoms': set(ring_atoms),
                    'bonds': set(bond_rings[ring_idx])
                })

    if not thf_rings:
        return False, "No tetrahydrofuran rings found"

    # For each tetrahydrofuran ring, check if it's fused to exactly two cyclopentane rings
    for thf_ring in thf_rings:
        fused_rings = []
        thf_bonds = thf_ring['bonds']

        # Find rings that share bonds with the THF ring (fused rings)
        for ring_idx, ring_bonds in enumerate(bond_rings):
            if ring_idx == thf_ring['idx']:
                continue  # Skip the THF ring itself
            # Check if the rings share any bonds (fused)
            shared_bonds = thf_bonds.intersection(ring_bonds)
            if shared_bonds:
                fused_rings.append({
                    'idx': ring_idx,
                    'atoms': set(atom_rings[ring_idx]),
                    'bonds': set(ring_bonds)
                })

        # Check if there are exactly two fused rings
        if len(fused_rings) != 2:
            continue  # Not the desired tricyclic system

        # Check if both fused rings are 5-membered rings composed only of carbon atoms
        valid_fused = True
        for fused_ring in fused_rings:
            # Check ring size
            if len(fused_ring['atoms']) != 5:
                valid_fused = False
                break
            # Check if all atoms are carbons
            if not all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in fused_ring['atoms']):
                valid_fused = False
                break

        if valid_fused:
            # Found a dihydroagarofuran skeleton
            return True, "Contains dihydroagarofuran skeleton with fused tricyclic system"

    # If no valid skeleton found
    return False, "Dihydroagarofuran skeleton not found"