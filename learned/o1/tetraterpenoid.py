"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26976 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is derived from a tetraterpene (C40 backbone), possibly modified
    by rearrangement or removal of methyl groups. Tetraterpenoids may also have
    common modifications like hydroxylation, epoxidation, glycosylation, or esterification.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple:
            bool: True if molecule is a tetraterpenoid, False otherwise
            str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Allowable elements for tetraterpenoids
    allowable_elements = {1, 6, 7, 8, 16, 15}  # H, C, N, O, S, P

    # Check for disallowed elements
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowable_elements:
            return False, f"Contains disallowed element: {atom.GetSymbol()}"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Adjusted carbon count range to include modifications (e.g., glycosides, cleavage)
    if c_count < 20 or c_count > 60:
        return False, f"Carbon count ({c_count}) not within typical range for tetraterpenoids"

    # Approximate number of isoprene units
    isoprene_units = c_count / 5
    if isoprene_units < 6:
        return False, f"Insufficient isoprene units (approximately {isoprene_units:.1f})"

    # Check for extended conjugation (conjugated double bonds)
    # Use RDKit's method to find the largest conjugated system
    largest_conj_len = 0
    for bond in mol.GetBonds():
        if bond.GetIsConjugated() and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Use a recursive method to find the length of the conjugated system
            visited_bonds = set()
            stack = [(bond, 1)]
            while stack:
                current_bond, length = stack.pop()
                if current_bond.GetIdx() in visited_bonds:
                    continue
                visited_bonds.add(current_bond.GetIdx())
                begin_atom = current_bond.GetBeginAtom()
                end_atom = current_bond.GetEndAtom()
                for neighbor_bond in begin_atom.GetBonds():
                    if neighbor_bond.GetIdx() != current_bond.GetIdx():
                        if neighbor_bond.GetIsConjugated() and neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            stack.append((neighbor_bond, length + 1))
                for neighbor_bond in end_atom.GetBonds():
                    if neighbor_bond.GetIdx() != current_bond.GetIdx():
                        if neighbor_bond.GetIsConjugated() and neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            stack.append((neighbor_bond, length + 1))
                if length > largest_conj_len:
                    largest_conj_len = length

    if largest_conj_len < 10:
        return False, f"Insufficient extended conjugation (maximum conjugated length is {largest_conj_len})"

    # Count number of oxygen atoms (allowing for common functional groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 20:
        return False, f"Too many oxygen atoms ({o_count}) for a typical tetraterpenoid"

    # Exclude molecules with uncommon functional groups
    disallowed_patterns = [
        Chem.MolFromSmarts("C(=O)N"),  # Amide bond
        Chem.MolFromSmarts("[!#1;!#6;!#7;!#8;!#15;!#16]"),  # Atoms other than H, C, N, O, P, S
        Chem.MolFromSmarts("[Mg,Ca,Fe,Cu,Zn]"),  # Metals
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups or elements"

    # Check molecular weight to ensure it's within a reasonable range
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1500:
        return False, f"Molecular weight ({mol_wt:.2f}) not typical for a tetraterpenoid"

    return True, "Molecule matches criteria for a tetraterpenoid (modified C40 backbone)"