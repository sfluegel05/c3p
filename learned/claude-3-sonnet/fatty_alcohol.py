"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: CHEBI:36195 fatty alcohol

A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than
27 carbon atoms. Fatty alcohols may be saturated or unsaturated and may be branched
or unbranched. They may contain one or more alcohol groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count alcohol groups (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) < 1:
        return False, "No alcohol groups found"

    # Find longest aliphatic carbon chain
    chain_pattern = Chem.MolFromSmarts("[CH2][CH2][CH2]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No aliphatic carbon chain found"

    chain_lengths = []
    for match in chain_matches:
        chain_start = match[0]
        chain = Chem.Mol.GetAtomBrandsDepthFirstSearch(mol, chain_start)
        chain_lengths.append(len(chain))

    longest_chain_length = max(chain_lengths)

    if longest_chain_length < 3 or longest_chain_length > 27:
        return False, f"Longest carbon chain length ({longest_chain_length}) outside of allowed range (3-27)"

    # Check for aliphatic (allow branched structures)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains cyclic structures, not allowed"

    # Count unsaturations (double bonds and triple bonds)
    num_double_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol) - rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(mol))
    num_triple_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol) - rdMolDescriptors.CalcNumRotatableBonds(Chem.RemoveHs(mol, implicitOnly=False))
    num_unsaturations = num_double_bonds + num_triple_bonds

    if num_unsaturations > longest_chain_length - 1:
        return False, "Too many unsaturations for carbon chain length"

    # Check for forbidden functional groups
    forbidden_patterns = ["[N]", "[S]", "[P]", "C(=O)O", "C(=O)N", "C(=O)C(=O)", "O=C-O", "O-C=O", "C-O-C"]
    for pattern in forbidden_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Contains forbidden functional group: {pattern}"

    return True, "Aliphatic alcohol with 3-27 carbon chain length and one or more alcohol groups"