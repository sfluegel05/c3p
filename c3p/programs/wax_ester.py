"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: CHEBI:50447 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of
    the carboxy group of a fatty acid with the alcoholic hydroxy group of
    a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester functional group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if len(ester_matches) == 0:
        return False, "No ester functional group found"
    elif len(ester_matches) > 1:
        return False, f"Multiple ester groups found ({len(ester_matches)}), molecule should have only one ester linkage"

    # Get the atoms involved in the ester group
    ester_atoms = ester_matches[0]
    carbonyl_c = ester_atoms[0]
    oxygen = ester_atoms[2]

    # Get the two chains attached to the ester
    # Chain from carbonyl carbon (fatty acid side)
    fatty_acid_chain = []
    visited = set()
    def traverse_chain(atom_idx):
        for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx != ester_atoms[1] and n_idx not in visited:
                visited.add(n_idx)
                if neighbor.GetSymbol() == 'C':
                    fatty_acid_chain.append(n_idx)
                    traverse_chain(n_idx)
    traverse_chain(carbonyl_c)

    fatty_acid_length = len(fatty_acid_chain)

    # Chain from alkoxy oxygen (fatty alcohol side)
    fatty_alcohol_chain = []
    visited = set()
    def traverse_chain_alcohol(atom_idx):
        for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx != ester_atoms[1] and n_idx not in visited:
                visited.add(n_idx)
                if neighbor.GetSymbol() == 'C':
                    fatty_alcohol_chain.append(n_idx)
                    traverse_chain_alcohol(n_idx)
    traverse_chain_alcohol(oxygen)

    fatty_alcohol_length = len(fatty_alcohol_chain)

    # Check if both chains are long enough to be considered fatty chains (typically ≥12 carbons)
    if fatty_acid_length < 12:
        return False, f"Fatty acid chain too short ({fatty_acid_length} carbons), need at least 12"
    if fatty_alcohol_length < 12:
        return False, f"Fatty alcohol chain too short ({fatty_alcohol_length} carbons), need at least 12"

    # Check for other functional groups (should be minimal in wax esters)
    num_func_groups = rdMolDescriptors.CalcNumFunctionalGroups(mol)
    num_rings = mol.GetRingInfo().NumRings()
    if num_func_groups > 1 or num_rings > 0:
        return False, "Molecule has additional functional groups or rings not characteristic of wax esters"

    return True, "Molecule is a wax ester with long-chain fatty acid and fatty alcohol"