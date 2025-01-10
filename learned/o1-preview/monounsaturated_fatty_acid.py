"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid (MUFA) is a fatty acid with exactly one double or triple bond
    in the fatty acid chain and singly bonded carbon atoms in the rest of the chain.
    MUFAs may contain rings or branching and can have additional functional groups like hydroxyls.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove salts and keep only the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, default=mol, key=lambda m: m.GetNumAtoms())

    # Identify carboxylic acid or ester group (allowing for esterified fatty acids)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1,-]')
    ester = Chem.MolFromSmarts('C(=O)O[CX4]')
    acid_or_ester = False

    if mol.HasSubstructMatch(carboxylic_acid):
        acid_or_ester = True
        is_ester = False
        fxn_group = carboxylic_acid
    elif mol.HasSubstructMatch(ester):
        acid_or_ester = True
        is_ester = True
        fxn_group = ester

    if not acid_or_ester:
        return False, "No carboxylic acid or ester group found"

    # Check total number of carboxylate or ester groups (should be one)
    num_acid_groups = len(mol.GetSubstructMatches(carboxylic_acid))
    num_ester_groups = len(mol.GetSubstructMatches(ester))
    total_acid_ester_groups = num_acid_groups + num_ester_groups
    if total_acid_ester_groups != 1:
        return False, f"Expected one carboxylic acid or ester group, found {total_acid_ester_groups}"

    # Count phosphate groups to exclude phospholipids
    phosphate = Chem.MolFromSmarts('P(=O)([OX1])[OX1]')
    if mol.HasSubstructMatch(phosphate):
        return False, "Molecule contains phosphate group, likely a phospholipid"

    # Check the molecule is not too large (exclude complex lipids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, "Molecular weight too high for a fatty acid"

    # Find the carbon chain(s) connected to the carboxyl carbon
    # Identify the carbonyl carbon atom
    matches = mol.GetSubstructMatches(fxn_group)
    carbonyl_c_idx = matches[0][0]  # first carbonyl carbon atom index
    carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

    # Use a breadth-first search to find the longest carbon chain starting from the carbonyl carbon
    from collections import deque

    visited = set()
    queue = deque()
    chain_atoms = set()
    double_bonds = 0
    triple_bonds = 0

    queue.append((carbonyl_c, None))
    while queue:
        atom, prev_atom = queue.popleft()
        atom_idx = atom.GetIdx()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)

        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_atoms.add(atom_idx)
        elif atom.GetAtomicNum() not in (1, 6, 8):  # Allow H, C, O atoms
            return False, "Molecule contains atoms other than carbon, hydrogen, and oxygen in the chain"

        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if prev_atom and nbr.GetIdx() == prev_atom.GetIdx():
                continue
            bond_type = bond.GetBondType()
            if bond_type == Chem.rdchem.BondType.DOUBLE:
                if atom.GetAtomicNum() == 6 and nbr.GetAtomicNum() == 6:
                    double_bonds += 1
            elif bond_type == Chem.rdchem.BondType.TRIPLE:
                if atom.GetAtomicNum() == 6 and nbr.GetAtomicNum() == 6:
                    triple_bonds += 1
            queue.append((nbr, atom))

    total_multiple_bonds = double_bonds + triple_bonds
    if total_multiple_bonds != 1:
        return False, f"Expected exactly one double or triple bond in the chain, found {total_multiple_bonds}"

    # Check for the number of carboxyl groups (should be one)
    carboxyl_groups = Chem.MolFromSmarts('C(=O)[O;H1,-]')
    num_carboxyl = len(mol.GetSubstructMatches(carboxyl_groups))
    if num_carboxyl != 1:
        return False, f"Expected one carboxyl group, found {num_carboxyl}"

    # Exclude molecules with multiple acyl chains (e.g., triglycerides)
    ester_oxygen = Chem.MolFromSmarts('[CX3](=O)O[CX3](=O)')
    if mol.HasSubstructMatch(ester_oxygen):
        return False, "Molecule contains multiple ester linkages, possible triglyceride or complex lipid"

    # Allow rings and branching in the chain, as some MUFAs have cyclic structures or branches
    # Ensure that all heavy atoms are connected (i.e., the molecule is contiguous)
    if not Chem.rdmolops.GetMolFrags(mol, sanitized=False, asMols=False, frags=...)[1] == 1:
        return False, "Molecule is not a single continuous structure"

    return True, "Molecule is a monounsaturated fatty acid"