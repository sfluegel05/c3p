"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:89998 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine consists of sphinganine with a fatty acyl group attached to the amino group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core sphinganine pattern with stereochemistry
    # Pattern matches: [NH]C(=O)-[C@...]-[C@...](O)-CH2OH and long carbon chain
    sphinganine_pattern = Chem.MolFromSmarts(
        '[NH1]C(=O)[C@H]([C@H](O)CCCCCCCCCCCCCCC)CO'  # Adjust chain length as needed
    )
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "Missing core sphinganine backbone with acyl group"

    # Verify fatty acyl chain length (at least 8 carbons in acyl group)
    amide_n = mol.GetSubstructMatch(Chem.MolFromSmarts('[NH1]C(=O)'))[0]
    acyl_start = mol.GetAtomWithIdx(amide_n).GetNeighbors()[1]  # Carbonyl carbon
    chain_length = 0
    visited = set()
    stack = [(acyl_start, 0)]
    
    while stack:
        atom, depth = stack.pop()
        if atom.GetIdx() in visited or atom.GetAtomicNum() != 6:
            continue
        visited.add(atom.GetIdx())
        chain_length += 1
        # Follow single bonds only, exclude branches
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6:
                    stack.append((neighbor, depth + 1))
    
    # Subtract 1 for the carbonyl carbon itself
    if (chain_length - 1) < 8:
        return False, f"Acyl chain too short (length {chain_length-1})"

    return True, "Contains N-acylsphinganine structure with proper acyl chain"