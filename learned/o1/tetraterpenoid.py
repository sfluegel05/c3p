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

    # Allowable heteroatoms for tetraterpenoids
    allowable_atoms = {1, 6, 7, 8, 15}  # H, C, N, O, P

    # Check for disallowed heteroatoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowable_atoms:
            return False, f"Contains heteroatom ({atom.GetSymbol()}) not typical for tetraterpenoids"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Adjusted carbon count range to include modifications (e.g., glycosides)
    if c_count < 35 or c_count > 60:
        return False, f"Carbon count ({c_count}) not within typical range for tetraterpenoids"

    # Identify isoprene units (C5 units)
    # Isoprene unit patterns allowing for modifications
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C-C")  # Simple isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)

    if num_isoprene_units < 6:
        return False, f"Insufficient isoprene units detected ({num_isoprene_units})"

    # Check for extended conjugation (conjugated double bonds)
    # Count the number of conjugated double bonds
    conjugated_bonds = 0
    bonds = mol.GetBonds()
    for bond in bonds:
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                # Check if connected to another double bond (conjugation)
                for neighbor in begin_atom.GetNeighbors():
                    if neighbor.GetIdx() != end_atom.GetIdx():
                        neighbor_bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), begin_atom.GetIdx())
                        if neighbor_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            conjugated_bonds += 1
                            break

    if conjugated_bonds < 10:
        return False, f"Insufficient extended conjugation (only {conjugated_bonds} conjugated bonds)"

    # Check for common tetraterpenoid functional groups
    # Allow hydroxyls, carbonyls, epoxides, glycosides, esters

    # Count number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # Allow oxygen count up to 20 to include glycosides, esters
    if o_count > 20:
        return False, f"Too many oxygen atoms ({o_count}) for a typical tetraterpenoid"

    # Exclude molecules with unusual functional groups
    # Disallowed functional groups (e.g., amino acids, peptides, metals)
    disallowed_patterns = [
        Chem.MolFromSmarts("C(=O)N"),  # Amide bond
        Chem.MolFromSmarts("[!#6;!#1]"),  # Any atom not carbon or hydrogen
    ]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed functional groups"

    # Calculate molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1500:
        return False, f"Molecular weight ({mol_wt:.2f}) not typical for a tetraterpenoid"

    return True, "Molecule matches criteria for a tetraterpenoid (C40 backbone with modifications)"