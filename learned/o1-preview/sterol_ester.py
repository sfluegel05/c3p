"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sterol nucleus pattern (cyclopentanoperhydrophenanthrene)
    sterol_nucleus_smarts = '[#6]12~[#6]3~[#6](~[#6]~[#6]~1)~[#6](~[#6]~[#6]~2)~[#6]~[#6]~[#6]3'
    sterol_nucleus = Chem.MolFromSmarts(sterol_nucleus_smarts)

    if not mol.HasSubstructMatch(sterol_nucleus):
        return False, "No sterol nucleus found"

    # Find all matches of sterol nucleus
    sterol_matches = mol.GetSubstructMatches(sterol_nucleus)
    if not sterol_matches:
        return False, "No sterol nucleus found"

    # For each match, check for esterified 3-hydroxy group
    ester_found = False
    for match in sterol_matches:
        sterol_atoms = set(match)
        mol_with_H = Chem.AddHs(mol)
        
        # Map the atoms in the match to the molecule with hydrogens
        atom_map = {}
        for idx, atom_idx in enumerate(match):
            atom_map[atom_idx] = mol_with_H.GetAtomWithIdx(atom_idx)

        # Assuming atom index 2 corresponds to the C3 position in sterols
        # Need to identify the 3-position hydroxyl group
        # Use a SMARTS pattern for 3-hydroxy sterol
        hydroxy_smarts = '[#6;r3]-[#6;r6](-[OX2H1])-[#6;r6]'
        hydroxy_pattern = Chem.MolFromSmarts(hydroxy_smarts)

        if mol.HasSubstructMatch(hydroxy_pattern):
            # Check if hydroxy group is esterified
            ester_smarts = '[#6]-[OX2H0]-C=O'
            ester_pattern = Chem.MolFromSmarts(ester_smarts)
            ester_bonds = mol.GetSubstructMatches(ester_pattern)

            # Check if ester oxygen is connected to the 3-position carbon
            for bond in ester_bonds:
                ester_oxygen = bond[1]  # Oxygen atom in ester group
                ester_carbon = bond[2]  # Carbonyl carbon in ester group

                # Check if the ester oxygen is connected to the 3-position carbon
                if ester_oxygen in sterol_atoms:
                    ester_found = True
                    break
            if ester_found:
                break
        else:
            continue

    if ester_found:
        return True, "Contains sterol nucleus with esterified 3-hydroxy group"
    else:
        return False, "No esterification at 3-hydroxy group of sterol nucleus found"