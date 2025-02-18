"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:15727 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a ribose/deoxyribose sugar with a phosphate group
    attached to the 5' carbon and a nitrogenous base (purine/pyrimidine) attached to the 1' carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group attached to 5' position of sugar
    # Pattern: Phosphate connected via ester to CH2 group attached to sugar ring oxygen
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)(O)(O)OCC1OC(CO1)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to 5' position of sugar"

    # Check for nitrogenous base (purine/pyrimidine) attached to sugar
    # Look for at least two nitrogen atoms in rings connected to the sugar
    base_nitrogens = 0
    rings = mol.GetRingInfo().AtomRings()
    for ring in rings:
        # Check if ring contains at least two nitrogens and is connected to the sugar
        ring_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if ring_nitrogens >= 2:
            # Check if any atom in this ring is connected to the sugar
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C' and neighbor.IsInRing():
                        # Check if the neighbor is part of the sugar ring
                        sugar_ring_atoms = [a for a in mol.GetAtomWithIdx(neighbor.GetIdx()).GetOwningMol().GetAtomWithIdx(neighbor.GetIdx()).GetNeighbors() if a.IsInRing()]
                        if any(a.GetSymbol() == 'O' for a in sugar_ring_atoms):
                            base_nitrogens += ring_nitrogens
                            break
                if base_nitrogens > 0:
                    break
            if base_nitrogens >= 2:
                break
    if base_nitrogens < 2:
        return False, "Insufficient nitrogen atoms in base rings"

    return True, "Contains 5'-phosphate group attached to ribose/deoxyribose with nitrogenous base"