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

    # Define the N-acylglycine pattern
    # Pattern represents R-C(=O)-N-CH2-C(=O)-O
    n_acylglycine_smarts = '[#6](=O)-[NX3]-[CH2]-C(=O)[O;H1,-]'
    n_acylglycine_pattern = Chem.MolFromSmarts(n_acylglycine_smarts)
    if n_acylglycine_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Attempt to find the N-acylglycine substructure
    matches = mol.GetSubstructMatches(n_acylglycine_pattern)
    if not matches:
        return False, "No N-acylglycine substructure found"

    # Iterate over matches to ensure correct classification
    for match in matches:
        # Atoms indices from the match
        # match = (acyl carbon, nitrogen, CH2 carbon, carboxyl carbon, oxygen)
        acyl_c_atom = mol.GetAtomWithIdx(match[0])
        n_atom = mol.GetAtomWithIdx(match[1])
        ch2_atom = mol.GetAtomWithIdx(match[2])
        carboxyl_c_atom = mol.GetAtomWithIdx(match[3])

        # Confirm that the CH2 carbon (glycine alpha carbon) has exactly two non-hydrogen neighbors
        ch2_neighbors = [atom for atom in ch2_atom.GetNeighbors() if atom.GetAtomicNum() != 1]
        if len(ch2_neighbors) != 2:
            continue  # Skip to next match

        # Confirm that the nitrogen is part of an amide linkage (connected to acyl carbon)
        if not acyl_c_atom.IsInRing():
            bond = mol.GetBondBetweenAtoms(acyl_c_atom.GetIdx(), n_atom.GetIdx())
            if bond is None or bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue  # Not an amide bond

        # Confirm that the carboxyl carbon is a carboxylic acid (connected to OH)
        o_atoms = [atom for atom in carboxyl_c_atom.GetNeighbors() if atom.GetAtomicNum() == 8]
        has_carboxylic_acid = False
        for o_atom in o_atoms:
            bond = mol.GetBondBetweenAtoms(carboxyl_c_atom.GetIdx(), o_atom.GetIdx())
            if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                if o_atom.GetTotalNumHs() == 1:
                    has_carboxylic_acid = True
                    break
        if not has_carboxylic_acid:
            continue  # Not a carboxylic acid

        # If all checks pass, molecule is an N-acylglycine
        return True, "Molecule is an N-acylglycine"

    return False, "Molecule does not match N-acylglycine criteria"