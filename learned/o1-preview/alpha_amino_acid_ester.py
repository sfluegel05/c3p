"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
"""

from rdkit import Chem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is the amino acid ester derivative obtained by the formal
    condensation of an alpha-amino acid with an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification or misclassification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    amine_pattern = Chem.MolFromSmarts("[NX3;H1,H2]")  # Primary or secondary amine
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CX4]")  # Ester group
    alpha_carbon_pattern = Chem.MolFromSmarts("[CX4H]")  # Tetrahedral carbon with at least one hydrogen

    # Find all amine nitrogens
    amine_atoms = mol.GetSubstructMatches(amine_pattern)
    if not amine_atoms:
        return False, "No primary or secondary amine found"

    # Find all ester carbons
    ester_atoms = mol.GetSubstructMatches(ester_pattern)
    if not ester_atoms:
        return False, "No ester group found"

    # Find all alpha carbons
    alpha_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' 
                     and atom.GetHybridization() == Chem.HybridizationType.SP3 
                     and atom.GetTotalNumHs() >= 1]

    # Check for alpha-amino acid ester connectivity
    for alpha_c_idx in alpha_carbons:
        alpha_c = mol.GetAtomWithIdx(alpha_c_idx)

        # Check if alpha carbon is connected to an amine nitrogen
        connected_to_amine = False
        for neighbor in alpha_c.GetNeighbors():
            if neighbor.GetIdx() in [idx for (idx,) in amine_atoms]:
                connected_to_amine = True
                break
        if not connected_to_amine:
            continue

        # Check if alpha carbon is connected to a carbonyl carbon of an ester group
        connected_to_ester = False
        for neighbor in alpha_c.GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                # Check if this carbon is part of an ester group
                for ester_carbon_idx, in ester_atoms:
                    if neighbor.GetIdx() == ester_carbon_idx:
                        connected_to_ester = True
                        break
            if connected_to_ester:
                break
        if connected_to_ester:
            return True, "Contains alpha-amino acid ester moiety"

    return False, "Does not contain alpha-amino acid ester moiety"


__metadata__ = {
    'chemical_class': {
        'name': 'alpha-amino acid ester',
        'definition': 'The amino acid ester derivative obtained the formal condensation of an alpha-amino acid with an alcohol.',
    },
}