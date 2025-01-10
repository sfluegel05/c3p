"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:57305 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is an amino acid where the carboxyl group is deprotonated to COO-
    and the amino group is protonated to NH3+; major species at physiological pH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for alpha-amino-acid zwitterion
    pattern = Chem.MolFromSmarts('[#6H1]-[#7+;H3]-[#6](=O)[O-]')
    if not mol.HasSubstructMatch(pattern):
        return False, "Missing alpha-amino-acid zwitterion motif"

    # Find all matches of the pattern
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) == 0:
        return False, "No alpha-amino-acid zwitterion center found"
    elif len(matches) > 1:
        return False, "Multiple alpha-amino-acid zwitterion centers found"

    # Verify that the molecule has no net charge
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != 0:
        return False, f"Molecule has net charge {total_charge}, should be neutral"

    # Check that the alpha carbon is connected to a side chain (any atom)
    match = matches[0]
    alpha_carbon_idx = match[0]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    neighbors = alpha_carbon.GetNeighbors()
    
    # Count the number of non-hydrogen, non-functional group neighbors (side chain)
    side_chain_atoms = [atom for atom in neighbors if atom.GetAtomicNum() != 1 and atom.GetIdx() not in match[1:]]
    if len(side_chain_atoms) == 0:
        return False, "Alpha carbon missing side chain"

    return True, "Molecule is an alpha-amino-acid zwitterion"

# Metadata
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:57305',
        'name': 'alpha-amino-acid zwitterion',
        'definition': 'An amino acid zwitterion obtained by transfer of a proton from the carboxy to the amino group of any alpha-amino acid; major species at pH 7.3.',
        'parents': ['CHEBI:33518', 'CHEBI:59703']
    }
}