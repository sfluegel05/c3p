"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide consists of an inositol residue linked via a phosphodiester bridge
    to a ceramide moiety (sphingoid base linked via an amide bond to a fatty acyl chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify inositol ring (six-membered ring with at least 5 hydroxyl groups attached)
    inositol_found = False
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) == 6:
            hydroxyl_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Check if atom is carbon
                if atom.GetAtomicNum() != 6:
                    break
                # Count oxygens attached to the carbon (excluding ring atoms)
                oxygen_attached = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and not neighbor.IsInRing():
                        oxygen_attached = True
                if oxygen_attached:
                    hydroxyl_count += 1
            # If at least 5 carbons in the ring have hydroxyl groups attached
            if hydroxyl_count >= 5:
                inositol_ring_idx = set(ring)
                inositol_found = True
                break
    if not inositol_found:
        return False, "No inositol ring found"

    # Identify phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Check connectivity between inositol ring and phosphate group
    phosphate_connected_to_inositol = False
    for match in phosphate_matches:
        phosphate_atom_idx = match[0]  # Index of phosphorus atom
        for bond in mol.GetAtomWithIdx(phosphate_atom_idx).GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(phosphate_atom_idx)
            if neighbor_idx in inositol_ring_idx:
                phosphate_connected_to_inositol = True
                phosphate_atom = mol.GetAtomWithIdx(phosphate_atom_idx)
                phosphate_neighbors = [n.GetIdx() for n in phosphate_atom.GetNeighbors()]
                break
        if phosphate_connected_to_inositol:
            break
    if not phosphate_connected_to_inositol:
        return False, "Phosphate group is not connected to inositol ring"

    # Identify ceramide moiety (amide bond linked to long-chain base)
    # Look for amide bond
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide bond found (ceramide moiety)"
    # Check if nitrogen is part of a long-chain sphingoid base
    ceramide_found = False
    for amide_match in amide_matches:
        nitrogen_idx = amide_match[2]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)
        # Get the chain connected to nitrogen (excluding the amide bond)
        sphingoid_chain = Chem.FragmentOnBonds(mol, [nitrogen_atom.GetBonds()[0].GetIdx()], addDummies=False)
        sphingoid_atom = sphingoid_chain.GetAtomWithIdx(nitrogen_idx)
        # Count carbon atoms connected in a chain (excluding rings)
        carbon_chain_length = 0
        visited = set()
        stack = [sphingoid_atom]
        while stack:
            current_atom = stack.pop()
            if current_atom.GetIdx() in visited:
                continue
            visited.add(current_atom.GetIdx())
            if current_atom.GetAtomicNum() == 6 and not current_atom.IsInRing():
                carbon_chain_length += 1
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() in [6, 8] and not neighbor.IsInRing():
                        stack.append(neighbor)
        if carbon_chain_length >= 12:
            ceramide_found = True
            break
    if not ceramide_found:
        return False, "No long-chain sphingoid base found (ceramide moiety)"

    # Check if phosphate group connects inositol and ceramide moieties
    phosphate_connections = phosphate_neighbors
    inositol_connected = False
    ceramide_connected = False
    for idx in phosphate_connections:
        if idx in inositol_ring_idx:
            inositol_connected = True
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 8:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and not neighbor.IsInRing():
                    # Possible connection to ceramide
                    ceramide_connected = True
    if not (inositol_connected and ceramide_connected):
        return False, "Phosphate group does not connect inositol and ceramide moieties"

    return True, "Contains inositol phosphoceramide structure"