"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:35641 polyprenol

A polyprenol is defined as: Any member of the class of prenols possessing the general formula
H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of more than one isoprene units.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal -OH group
    has_terminal_oh = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    if not has_terminal_oh:
        return False, "No terminal -OH group found"

    # Find the longest carbon chain
    longest_chain = get_longest_carbon_chain(mol)
    if longest_chain is None:
        return False, "No suitable carbon chain found"

    # Check if the carbon chain follows the polyprenol pattern
    is_valid_chain, reason = is_valid_polyprenol_chain(longest_chain)
    if not is_valid_chain:
        return False, reason

    # Optionally, check molecular weight range or other properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for polyprenols"

    return True, "Molecule contains a carbon chain following the polyprenol pattern, with a terminal -OH group"

def get_longest_carbon_chain(mol):
    """
    Returns the longest carbon chain in the molecule.
    """
    longest_chain = None
    max_length = 0

    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() != 6 or atom2.GetAtomicNum() != 6:
            continue

        chain = [atom1, atom2]
        while True:
            neighbors = [neighbor for neighbor in chain[-1].GetNeighbors() if neighbor.GetAtomicNum() == 6 and neighbor not in chain]
            if not neighbors:
                break
            chain.append(neighbors[0])

        if len(chain) > max_length:
            longest_chain = chain
            max_length = len(chain)

    return longest_chain

def is_valid_polyprenol_chain(chain):
    """
    Checks if the given carbon chain follows the polyprenol pattern: C-C=C-C-C=C-C-...
    with at least two instances of the C=C-C pattern and the correct number of methyl groups.
    """
    pattern_count = 0
    methyl_count = 0

    for i in range(len(chain) - 2):
        atom1, atom2, atom3 = chain[i:i+3]
        if atom1.GetHybridization() == Chem.HybridizationType.SP2 and atom2.GetHybridization() == Chem.HybridizationType.SP3 and atom3.GetHybridization() == Chem.HybridizationType.SP2:
            pattern_count += 1
            methyl_count += sum(1 for neighbor in atom2.GetNeighbors() if neighbor.GetAtomicNum() == 6 and neighbor.GetTotalNumHs() == 3)

    if pattern_count < 2:
        return False, "Less than two isoprene units found in the carbon chain"
    if methyl_count != pattern_count:
        return False, "Incorrect number of methyl groups attached to the carbon chain"

    return True, "Carbon chain follows the polyprenol pattern"