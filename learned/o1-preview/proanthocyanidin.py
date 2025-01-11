"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavan-3-ol substructure pattern
    # This SMARTS pattern represents the core flavan-3-ol structure
    flavan3ol_smarts = """
    [#6]1(:[#6]:[#6]:[#6]:[#6]:[#6]:1)              # Aromatic A ring
    [C@@H]2                                         # Chiral center at C2
    ([C@H](O)C3=CC=CC=C3)                           # Heterocyclic C ring with hydroxyl at C3, attached to aromatic B ring
    O[C@H]2                                         # Oxygen in ring closure
    """
    flavan3ol_smarts = flavan3ol_smarts.replace('\n', '').replace('    ', '')
    flavan3ol_pattern = Chem.MolFromSmarts(flavan3ol_smarts)
    if flavan3ol_pattern is None:
        return False, "Invalid flavan-3-ol SMARTS pattern"

    # Find flavan-3-ol units in the molecule
    matches = mol.GetSubstructMatches(flavan3ol_pattern, useChirality=True)
    num_units = len(matches)

    if num_units < 2:
        return False, f"Contains only {num_units} flavan-3-ol unit(s), need at least 2"

    # Map each match to the atom indices corresponding to key positions
    # For flavan-3-ol units, we can attempt to identify C4, C6, and C8 positions for linkage
    # However, without explicit labels, this is complex
    # Instead, we check for bonds between units

    # Collect sets of atom indices for each flavan-3-ol unit
    unit_atom_sets = [set(match) for match in matches]

    # Check for linkages between units
    linked_units = 0
    for i in range(len(unit_atom_sets)):
        for j in range(i+1, len(unit_atom_sets)):
            # Check if there's a bond between any atom of unit i and unit j
            connected = False
            for atom_idx_i in unit_atom_sets[i]:
                for atom_idx_j in unit_atom_sets[j]:
                    bond = mol.GetBondBetweenAtoms(atom_idx_i, atom_idx_j)
                    if bond is not None:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                linked_units += 1

    if linked_units < 1:
        return False, "Flavan-3-ol units are not connected via typical linkages"

    # Additional checks can include molecular weight and rotatable bonds
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a proanthocyanidin"

    return True, f"Contains {num_units} flavan-3-ol units connected via typical linkages"