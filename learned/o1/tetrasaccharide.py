"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:28053 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    ring_atoms = ring_info.AtomRings()

    monosaccharide_units = []

    for ring in ring_atoms:
        ring_size = len(ring)
        if ring_size not in [5, 6]:
            continue  # Only consider furanose or pyranose rings

        # Collect atoms in ring
        ring_atom_indices = set(ring)
        ring_atoms_list = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Check if ring contains exactly one oxygen atom
        o_in_ring = [atom for atom in ring_atoms_list if atom.GetAtomicNum() == 8]
        if len(o_in_ring) != 1:
            continue

        # Check that all ring atoms are either carbon or oxygen
        if any(atom.GetAtomicNum() not in [6, 8] for atom in ring_atoms_list):
            continue

        # Identify ring oxygen and ring carbons
        ring_oxygen_idx = o_in_ring[0].GetIdx()
        ring_carbon_idxs = [atom.GetIdx() for atom in ring_atoms_list if atom.GetAtomicNum() == 6]

        # Check exocyclic oxygens on ring carbons
        exocyclic_oxygen_counts = 0
        possible_anomeric_carbons = []
        for idx in ring_carbon_idxs:
            atom = mol.GetAtomWithIdx(idx)
            oxygen_neighbors = 0
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == ring_oxygen_idx:
                    continue  # Skip the ring oxygen
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(idx, nbr_idx)
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        oxygen_neighbors +=1
            if oxygen_neighbors == 1:
                exocyclic_oxygen_counts +=1
            elif oxygen_neighbors == 2:
                # Anomeric carbon connected to two exocyclic oxygens
                possible_anomeric_carbons.append(idx)
            else:
                # Ring carbon without exocyclic oxygen, not a monosaccharide
                break
        else:
            # Ring passes checks
            monosaccharide_units.append({
                'ring_atoms': ring_atom_indices,
                'anomeric_carbons': possible_anomeric_carbons
            })
    
    # Remove duplicate rings (in case of shared rings)
    unique_rings = []
    for unit in monosaccharide_units:
        if unit['ring_atoms'] not in [u['ring_atoms'] for u in unique_rings]:
            unique_rings.append(unit)
    monosaccharide_units = unique_rings

    # Check if there are exactly four monosaccharide units
    if len(monosaccharide_units) != 4:
        return False, f"Found {len(monosaccharide_units)} monosaccharide units, need exactly 4"

    # Check for glycosidic linkages between the monosaccharide units
    glycosidic_bonds = 0
    for i, unit1 in enumerate(monosaccharide_units):
        for j, unit2 in enumerate(monosaccharide_units):
            if i >= j:
                continue  # Avoid double counting and self-comparison
            # Check for bonds between anomeric carbon of unit1 and exocyclic oxygen of unit2
            for anomeric_idx in unit1['anomeric_carbons']:
                anomeric_atom = mol.GetAtomWithIdx(anomeric_idx)
                for nbr in anomeric_atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx in unit2['ring_atoms']:
                        bond = mol.GetBondBetweenAtoms(anomeric_idx, nbr_idx)
                        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            # Found bond between anomeric carbon of unit1 and ring atom of unit2
                            glycosidic_bonds +=1
                            break

    # For a tetrasaccharide, there should be at least 3 glycosidic bonds connecting the 4 units
    if glycosidic_bonds < 3:
        return False, f"Found {glycosidic_bonds} glycosidic bonds, need at least 3"

    return True, "Contains four monosaccharide units linked via glycosidic bonds"