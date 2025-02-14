"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:18154 polyamine
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    Amino groups include primary, secondary, and tertiary amines (both aliphatic and aromatic),
    excluding nitrogen atoms in amides, nitro groups, nitriles, azo compounds, and other non-amino functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for amino groups
    # Primary amine: Nitrogen with two hydrogens and single bonds
    primary_amine = Chem.MolFromSmarts("[NX3;H2][#6]")
    # Secondary amine: Nitrogen with one hydrogen and single bonds
    secondary_amine = Chem.MolFromSmarts("[NX3;H1]([#6])[#6]")
    # Tertiary amine: Nitrogen with no hydrogens and single bonds
    tertiary_amine = Chem.MolFromSmarts("[NX3;H0]([#6])[#6]")
    # Aromatic amine (aniline type): Nitrogen attached to aromatic ring
    aromatic_amine = Chem.MolFromSmarts("[NX3;H1,H0][a]")

    # Combine all amine patterns
    amine_patterns = [primary_amine, secondary_amine, tertiary_amine, aromatic_amine]

    # Initialize set to keep track of unique amino nitrogen atoms
    amino_nitrogens = set()

    # Search for amino groups
    for pattern in amine_patterns:
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            n_idx = match[0]
            atom = mol.GetAtomWithIdx(n_idx)

            # Exclude nitrogens in amides (N-C=O)
            is_amide = False
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 6:  # Carbon
                    # Check if carbon is double bonded to oxygen
                    for nbr_bond in nbr.GetBonds():
                        nbr_atom = nbr_bond.GetOtherAtom(nbr)
                        if (nbr_bond.GetBondType() == rdchem.BondType.DOUBLE and
                            nbr_atom.GetAtomicNum() == 8):  # Oxygen
                            is_amide = True
                            break
                    if is_amide:
                        break
            if is_amide:
                continue

            # Exclude nitro groups (N(=O)-O)
            is_nitro = False
            oxy_count = 0
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:  # Oxygen
                    oxy_count += 1
            if oxy_count >= 2:
                is_nitro = True
            if is_nitro:
                continue

            # Exclude nitriles (N#C)
            is_nitrile = False
            for bond in atom.GetBonds():
                if bond.GetBondType() == rdchem.BondType.TRIPLE:
                    is_nitrile = True
                    break
            if is_nitrile:
                continue

            # Exclude azo compounds (N=N)
            is_azo = False
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if (nbr.GetAtomicNum() == 7 and
                    bond.GetBondType() == rdchem.BondType.DOUBLE):
                    is_azo = True
                    break
            if is_azo:
                continue

            # Passed all checks, add to amino nitrogens set
            amino_nitrogens.add(n_idx)

    amino_count = len(amino_nitrogens)

    if amino_count >= 2:
        return True, f"Molecule contains {amino_count} amino groups"
    else:
        return False, f"Molecule contains {amino_count} amino group(s), need at least 2"