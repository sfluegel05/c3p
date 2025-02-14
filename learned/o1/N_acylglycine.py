"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is an N-acyl-amino acid in which the amino acid specified is glycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the N-acylglycine SMARTS pattern
    # Pattern represents R-C(=O)-N-CH2-C(=O)-O
    n_acylglycine_smarts = '[#6;X3](=O)-N-CH2-C(=O)-O'
    n_acylglycine_pattern = Chem.MolFromSmarts(n_acylglycine_smarts)
    if n_acylglycine_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Attempt to find the N-acylglycine substructure
    matches = mol.GetSubstructMatches(n_acylglycine_pattern)
    if not matches:
        return False, "No N-acylglycine substructure found"

    # Verify each match to ensure correct classification
    for match in matches:
        # Extract atoms from the match
        acyl_c_atom = mol.GetAtomWithIdx(match[0])
        n_atom = mol.GetAtomWithIdx(match[1])
        ch2_atom = mol.GetAtomWithIdx(match[2])
        carboxyl_c_atom = mol.GetAtomWithIdx(match[3])
        o_atom = mol.GetAtomWithIdx(match[4])

        # Check that the nitrogen is part of an amide linkage
        bonds_to_n = [bond for bond in n_atom.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.SINGLE]
        if len(bonds_to_n) < 2:
            continue  # Nitrogen should be connected to acyl carbon and CH2

        # Confirm that the CH2 group is connected to both nitrogen and carboxyl carbon
        ch2_neighbors = [atom.GetIdx() for atom in ch2_atom.GetNeighbors()]
        if n_atom.GetIdx() not in ch2_neighbors or carboxyl_c_atom.GetIdx() not in ch2_neighbors:
            continue

        # Confirm that the carboxyl carbon is double-bonded to an oxygen
        carboxyl_bonds = carboxyl_c_atom.GetBonds()
        double_bonded_oxygens = [
            bond for bond in carboxyl_bonds
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetOtherAtom(carboxyl_c_atom).GetAtomicNum() == 8
        ]
        if not double_bonded_oxygens:
            continue

        # Confirm that the carboxyl carbon is single-bonded to a hydroxyl oxygen
        single_bonded_oxygen = [
            bond for bond in carboxyl_bonds
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and
            bond.GetOtherAtom(carboxyl_c_atom).GetAtomicNum() == 8 and
            bond.GetOtherAtom(carboxyl_c_atom).GetTotalDegree() == 1
        ]
        if not single_bonded_oxygen:
            continue

        # If all checks pass, molecule is an N-acylglycine
        return True, "Molecule is an N-acylglycine"

    return False, "Molecule does not match N-acylglycine criteria"