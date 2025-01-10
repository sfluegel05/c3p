"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify purine or pyrimidine base (allowing for substitutions)
    # Purine core: two fused rings with four nitrogens
    purine_base = Chem.MolFromSmarts('c1nc[nH0,n]c2ncnc12')
    # Pyrimidine core: six-membered ring with two nitrogens
    pyrimidine_base = Chem.MolFromSmarts('c1cncnc1')

    has_purine = mol.HasSubstructMatch(purine_base)
    has_pyrimidine = mol.HasSubstructMatch(pyrimidine_base)

    if not (has_purine or has_pyrimidine):
        return False, "Nucleobase not found (purine or pyrimidine base)"

    # Identify sugar ring (ribose or deoxyribose with possible modifications)
    # Five-membered ring with an oxygen and four carbons
    sugar = Chem.MolFromSmarts('[C;R1]1-[C;R1]-[O;R1]-[C;R1]-[C;R1]1')
    sugar_matches = mol.GetSubstructMatches(sugar)
    if not sugar_matches:
        return False, "Sugar ring (ribose or deoxyribose) not found"

    # Find the anomeric carbon (connected to both the sugar oxygen and the nucleobase)
    anomeric_carbons = []
    for match in sugar_matches:
        ring_atoms = set(match)
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
                if any(mol.GetAtomWithIdx(nbr_idx).GetAtomicNum() == 7 for nbr_idx in neighbors):
                    # Atom is bonded to nitrogen (possible glycosidic bond)
                    if any(nbr_idx in ring_atoms for nbr_idx in neighbors if mol.GetAtomWithIdx(nbr_idx).GetAtomicNum() == 8):
                        # Atom is also bonded to oxygen in ring (anomeric carbon)
                        anomeric_carbons.append(atom_idx)
    if not anomeric_carbons:
        return False, "N-glycosidic bond between nucleobase and sugar not found"

    # Check for phosphate group(s) attached to the 5' carbon of the sugar
    # Phosphate group pattern (allowing for multiple phosphates)
    phosphate = Chem.MolFromSmarts('OP(=O)(O)O')
    phosphate_matches = mol.GetSubstructMatches(phosphate)
    if not phosphate_matches:
        return False, "Phosphate group not found"

    # Identify the 5' carbon: carbon attached to sugar ring and to phosphate group
    phosphate_bonded = False
    for match in sugar_matches:
        ring_atoms = set(match)
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:
                # Look for carbon bonded to oxygen outside the ring (possible 5' carbon)
                for nbr in atom.GetNeighbors():
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in ring_atoms:
                        if nbr.GetAtomicNum() == 8:
                            # Check if oxygen is connected to phosphate
                            for nbr2 in nbr.GetNeighbors():
                                if nbr2.GetAtomicNum() == 15:
                                    phosphate_bonded = True
                                    break
                    if phosphate_bonded:
                        break
            if phosphate_bonded:
                break
        if phosphate_bonded:
            break

    if not phosphate_bonded:
        return False, "Phosphate group at 5' position not found"

    # Count the number of phosphate groups attached (mono-, di-, tri-, or tetra-phosphate)
    phosphorous_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    num_phosphates = len(phosphorous_atoms)
    if num_phosphates < 1 or num_phosphates > 4:
        return False, f"Found {num_phosphates} phosphate groups, but need between 1 and 4"

    return True, "Molecule is a nucleoside 5'-phosphate"

__metadata__ = {
    'chemical_class': {
        'name': "nucleoside 5'-phosphate",
        'definition': "A ribosyl or deoxyribosyl derivative of a pyrimidine or purine base in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated."
    }
}