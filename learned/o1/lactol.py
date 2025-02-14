"""
Classifies: CHEBI:38131 lactol
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by the intramolecular addition of a hydroxy group
    to an aldehydic or ketonic carbonyl group, resulting in a 1-oxacycloalkan-2-ol or unsaturated analogue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found_lactol = False

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Iterate over rings
    for ring in atom_rings:
        # Check if ring contains an oxygen atom
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_oxygen_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]

        for o_atom in ring_oxygen_atoms:
            # Get neighboring carbons in the ring
            for nbr in o_atom.GetNeighbors():
                if nbr.GetIdx() in ring and nbr.GetAtomicNum() == 6:
                    # Check if carbon is saturated (sp3 hybridized)
                    if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                        continue

                    # Check if carbon is bonded to exactly one oxygen (the ring oxygen)
                    o_count = sum(1 for atom in nbr.GetNeighbors() if atom.GetAtomicNum() == 8)
                    if o_count != 1:
                        continue

                    # Check if carbon has a hydroxyl group attached
                    has_oh = False
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 8 and nbr2.GetIdx() != o_atom.GetIdx():
                            # Check if this oxygen is a hydroxyl (bonded to hydrogen)
                            if nbr2.GetTotalHydrogenCount() == 1 and nbr2.IsInRing() == False:
                                has_oh = True
                            else:
                                has_oh = False
                                break  # Exclude if bonded to other oxygens or in ring
                    if has_oh:
                        # Exclude sugars by checking for multiple hydroxyl groups on carbons in the ring
                        ring_oxygen_counts = []
                        for atom in ring_atoms:
                            if atom.GetAtomicNum() == 6:
                                oh_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.IsInRing() == False)
                                ring_oxygen_counts.append(oh_count)
                        # If more than two carbons in the ring have hydroxyl groups, it's likely a sugar
                        if sum(1 for count in ring_oxygen_counts if count >=1) > 2:
                            continue

                        found_lactol = True
                        return True, "Contains lactol functional group"

    if not found_lactol:
        return False, "Lactol functional group not found"