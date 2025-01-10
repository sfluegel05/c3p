"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

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

    # Define the SMARTS pattern for O-acyl-L-carnitine
    # This pattern matches:
    # - A chiral carbon with S configuration [C@H]
    # - Connected to:
    #   - An esterified oxygen [O][C](=O)[*]
    #   - A carbon chain leading to a carboxylate group [CH2][C](=O)[O-]
    #   - A carbon chain leading to a quaternary ammonium group [CH2][N+](C)(C)C
    smarts = "[C@H](OC(=O)[#6])C[C](=O)[O-]"

    # Create the SMARTS molecule
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)

    # If no matches found, return False
    if not matches:
        return False, "O-acyl-L-carnitine pattern not found"

    # For each match, check the quaternary ammonium group
    for match in matches:
        chiral_atom_idx = match[0]
        chiral_atom = mol.GetAtomWithIdx(chiral_atom_idx)

        # Verify that the chiral center has S configuration
        if chiral_atom.HasProp('_CIPCode') and chiral_atom.GetProp('_CIPCode') != 'S':
            continue  # Not S configuration

        # Check for the quaternary ammonium group connected via two carbons
        connected = False
        for nbr in chiral_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Carbon
                # Look for CH2 connected to chiral center
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == chiral_atom_idx:
                        continue
                    if nbr2.GetAtomicNum() == 6:
                        # Look for carbon connected to nitrogen
                        for nbr3 in nbr2.GetNeighbors():
                            if nbr3.GetAtomicNum() == 7 and nbr3.GetFormalCharge() == 1:
                                # Check if nitrogen is quaternary ammonium [N+](C)(C)C
                                if len(nbr3.GetNeighbors()) == 4:
                                    methyl_count = sum(1 for atom in nbr3.GetNeighbors() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 1)
                                    if methyl_count == 3:
                                        connected = True
                                        break
                        if connected:
                            break
                if connected:
                    break
        if not connected:
            continue  # Quaternary ammonium group not found

        # All checks passed
        return True, "Molecule is an O-acyl-L-carnitine"

    # If none of the matches passed all checks
    return False, "Molecule does not match all criteria for O-acyl-L-carnitine"