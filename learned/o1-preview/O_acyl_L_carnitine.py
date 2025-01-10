"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is an O-acylcarnitine in which the carnitine component has L-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define SMARTS patterns for key features
    # Chiral carbon with L-configuration (S) connected to:
    # 1. Esterified oxygen with acyl group
    # 2. Quaternary ammonium group [N+](C)(C)C
    # 3. Carboxylate group [C](=O)[O-]

    # Pattern for chiral center with attached groups
    # Using [C@H] for L-configuration
    chiral_center_pattern = "[C@H]"

    # Ester linkage at the hydroxyl group (O-acyl)
    ester_pattern = "O[C](=O)[#6]"  # Oxygen connected to C=O and carbon chain

    # Quaternary ammonium group connected via two carbons
    quaternary_ammonium_pattern = "C[C+](C)(C)C"  # Trimethylammonium group

    # Carboxylate group connected via one carbon
    carboxylate_pattern = "C(=O)[O-]"

    # Combine patterns to match the backbone of O-acyl-L-carnitine
    # Chiral center connected to ester oxygen, N+ group, and carboxylate
    o_acyl_l_carnitine_smarts = f"{chiral_center_pattern}({ester_pattern})([#6]{quaternary_ammonium_pattern})[#6]{carboxylate_pattern}"

    # Create the SMARTS molecule
    pattern = Chem.MolFromSmarts(o_acyl_l_carnitine_smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)

    if not matches:
        return False, "O-acyl-L-carnitine pattern not found"

    for match in matches:
        chiral_atom_idx = match[0]
        chiral_atom = mol.GetAtomWithIdx(chiral_atom_idx)

        # Verify that the chiral center has S configuration
        if chiral_atom.HasProp('_CIPCode') and chiral_atom.GetProp('_CIPCode') != 'S':
            continue  # Not S configuration

        # Verify connections to ester oxygen, quaternary ammonium, and carboxylate groups
        ester_connected = False
        ammonium_connected = False
        carboxylate_connected = False

        for neighbor in chiral_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(chiral_atom_idx, neighbor_idx)
            bond_type = bond.GetBondType()

            if neighbor.GetAtomicNum() == 8:  # Oxygen (possible ester linkage)
                # Check if oxygen is part of ester linkage
                for ester_match in mol.GetSubstructMatches(Chem.MolFromSmarts(ester_pattern)):
                    if neighbor_idx in ester_match and chiral_atom_idx in ester_match:
                        ester_connected = True
                        break

            elif neighbor.GetAtomicNum() == 6:  # Carbon
                # Check for quaternary ammonium group
                for path in Chem.FindAllPathsOfLengthN(mol, 2, useBonds=False, useHs=False):
                    if chiral_atom_idx == path[0] and neighbor_idx == path[1]:
                        intermediate_atom = mol.GetAtomWithIdx(path[1])
                        if intermediate_atom.GetAtomicNum() == 6:
                            for nbr in intermediate_atom.GetNeighbors():
                                if nbr.GetAtomicNum() == 7 and nbr.GetFormalCharge() == 1:
                                    # Check if nitrogen is quaternary ammonium
                                    attached_carbons = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() == 6]
                                    if len(attached_carbons) == 3:
                                        ammonium_connected = True
                                        break
                            if ammonium_connected:
                                break

                # Check for carboxylate group
                for carboxylate_match in mol.GetSubstructMatches(Chem.MolFromSmarts(carboxylate_pattern)):
                    if neighbor_idx in carboxylate_match:
                        carboxylate_connected = True
                        break

        if ester_connected and ammonium_connected and carboxylate_connected:
            return True, "Molecule is an O-acyl-L-carnitine"

    return False, "Molecule does not match all criteria for O-acyl-L-carnitine"