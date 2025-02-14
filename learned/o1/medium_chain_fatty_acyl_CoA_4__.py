"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA molecule with an aliphatic acyl chain of 6-12 carbons,
    containing no rings, no aromatic systems, and composed only of C, H, and O atoms in the acyl chain.
    The acyl chain should not contain additional functional groups like carboxylic acids, phosphates, or heteroatoms.
    The phosphate and diphosphate groups are deprotonated (negative charges).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update property cache for the molecule
    mol.UpdatePropertyCache()

    # Check for phosphorus atoms (CoA derivative)
    contains_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
    if not contains_phosphorus:
        return False, "No phosphorus atoms found; molecule is not a CoA derivative"

    # Check for adenine moiety (part of CoA)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2N')
    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    if not adenine_matches:
        return False, "No adenine moiety found; molecule is not a CoA derivative"

    # Look for thioester linkage C(=O)-S-
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Use the first thioester match
    acyl_carbon_idx = thioester_matches[0][0]  # Carbonyl carbon
    sulfur_idx = thioester_matches[0][2]       # Sulfur atom

    # Collect acyl chain atoms starting from the carbonyl carbon, excluding sulfur atom
    acyl_chain_atoms = set()
    visited_atoms = set([sulfur_idx])  # Exclude sulfur atom

    current_atom = mol.GetAtomWithIdx(acyl_carbon_idx)
    neighbors = current_atom.GetNeighbors()

    # Exclude sulfur and oxygen atoms
    excluded_atom_idxs = [nbr.GetIdx() for nbr in neighbors if nbr.GetAtomicNum() in [8, 16]]  # Oxygen and Sulfur
    visited_atoms.update(excluded_atom_idxs)
    visited_atoms.add(acyl_carbon_idx)  # Exclude carbonyl carbon

    # Get the neighbor that is not sulfur or oxygen (start of acyl chain)
    acyl_chain_neighbors = [nbr for nbr in neighbors if nbr.GetIdx() not in visited_atoms]

    if len(acyl_chain_neighbors) != 1:
        return False, "Could not find acyl chain attached to carbonyl carbon"

    acyl_chain_start_idx = acyl_chain_neighbors[0].GetIdx()

    # Now, collect acyl chain atoms
    atoms_to_visit = [acyl_chain_start_idx]
    while atoms_to_visit:
        current_atom_idx = atoms_to_visit.pop()
        if current_atom_idx in visited_atoms:
            continue
        visited_atoms.add(current_atom_idx)
        acyl_chain_atoms.add(current_atom_idx)
        current_atom = mol.GetAtomWithIdx(current_atom_idx)
        for neighbor in current_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited_atoms:
                atoms_to_visit.append(neighbor_idx)

    # Create acyl chain fragment
    acyl_chain_frag = Chem.PathToSubmol(mol, list(acyl_chain_atoms))

    if acyl_chain_frag is None or acyl_chain_frag.GetNumAtoms() == 0:
        return False, "Could not isolate acyl chain fragment"

    # Update property cache for the fragment to initialize ring info
    acyl_chain_frag.UpdatePropertyCache()

    # Check for rings in acyl chain
    ring_info = acyl_chain_frag.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Acyl chain contains rings, not a fatty acyl chain"

    # Check for aromatic atoms in acyl chain
    if any(atom.GetIsAromatic() for atom in acyl_chain_frag.GetAtoms()):
        return False, "Acyl chain contains aromatic atoms, not a fatty acyl chain"

    # Check for allowed atoms (C, H, O) in acyl chain
    allowed_atoms = {6, 1, 8}  # Carbon, Hydrogen, Oxygen
    disallowed_atoms = set()
    for atom in acyl_chain_frag.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in allowed_atoms:
            disallowed_atoms.add(atom.GetSymbol())
    if disallowed_atoms:
        return False, f"Acyl chain contains disallowed atoms: {', '.join(disallowed_atoms)}"

    # Check for additional functional groups in acyl chain (carboxylic acids, phosphates)
    functional_group_patterns = [
        Chem.MolFromSmarts('C(=O)[O-]'),   # Carboxylate
        Chem.MolFromSmarts('P(=O)([O-])[O-]'),  # Phosphate
        Chem.MolFromSmarts('S'),  # Sulfur
        Chem.MolFromSmarts('N')   # Nitrogen
    ]
    for pattern in functional_group_patterns:
        if acyl_chain_frag.HasSubstructMatch(pattern):
            return False, "Acyl chain contains disallowed functional groups"

    # Count carbons in acyl chain
    carbon_count = sum(1 for atom in acyl_chain_frag.GetAtoms() if atom.GetAtomicNum() == 6)

    if 6 <= carbon_count <= 12:
        return True, f"Acyl chain has {carbon_count} carbons and is aliphatic; within medium-chain length (6-12)"
    else:
        return False, f"Acyl chain has {carbon_count} carbons, not within medium-chain length (6-12)"