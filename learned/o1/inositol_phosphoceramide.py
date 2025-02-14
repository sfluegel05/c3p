"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide consists of an inositol residue (possibly substituted) linked via
    a phosphodiester bridge to a ceramide moiety (sphingoid base linked via an amide bond to a fatty acyl chain).

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

    # Identify phosphate groups (phosphoric acid esters with two singly bonded oxygens)
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"

    # Identify inositol ring pattern (cyclohexane ring with at least 5 hydroxyl groups)
    inositol_pattern = Chem.MolFromSmarts("C1C[C@H](O)C[C@H](O)C[C@H](O)C1")
    # Allow for substitutions on the ring (e.g., mannose attached)
    # We'll check for a cyclohexane ring with at least 5 hydroxyl groups

    # Identify ceramide pattern (amide bond connected to a long chain with hydroxyl groups)
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H]([*])[#6]-[#6]-[#6]-[#6]-[#6]-[#6]")  # Simplified pattern

    for match in phosphate_matches:
        phosphorus_idx = match[0]
        phosphorus_atom = mol.GetAtomWithIdx(phosphorus_idx)

        # Get oxygen atoms singly bonded to phosphorus
        oxygen_atoms = [nbr for nbr in phosphorus_atom.GetNeighbors()
                        if nbr.GetAtomicNum() == 8 and
                        mol.GetBondBetweenAtoms(phosphorus_idx, nbr.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE]

        if len(oxygen_atoms) < 2:
            continue  # Not a phosphodiester bridge

        connected_fragments = []
        for oxygen_atom in oxygen_atoms:
            frag_atoms = rdmolops.GetConnectedAtomIndices(mol, oxygen_atom.GetIdx(), excludeAtoms={phosphorus_idx})
            frag_mol = Chem.PathToSubmol(mol, frag_atoms)
            connected_fragments.append(frag_mol)

        if len(connected_fragments) < 2:
            continue

        inositol_found = False
        ceramide_found = False

        for frag in connected_fragments:
            # Check for inositol ring with at least 5 hydroxyl groups
            ring_info = frag.GetRingInfo()
            has_inositol_ring = False
            for ring in ring_info.AtomRings():
                if len(ring) == 6:
                    # Count hydroxyl groups attached to ring carbons
                    hydroxyl_count = 0
                    for idx in ring:
                        atom = frag.GetAtomWithIdx(idx)
                        if atom.GetAtomicNum() == 6:  # Carbon
                            for nbr in atom.GetNeighbors():
                                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                                    hydroxyl_count += 1
                    if hydroxyl_count >= 5:
                        has_inositol_ring = True
                        break
            if has_inositol_ring:
                inositol_found = True
                continue

            # Check for ceramide moiety
            if frag.HasSubstructMatch(ceramide_pattern):
                ceramide_found = True
                continue

        if inositol_found and ceramide_found:
            return True, "Contains inositol phosphoceramide structure"

    return False, "Does not contain inositol phosphoceramide structure"