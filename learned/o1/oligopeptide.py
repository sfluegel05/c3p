"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids,
    typically between 2 and 20 residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond SMARTS pattern
    # This pattern specifically matches peptide bonds in the backbone
    peptide_bond_smarts = '[C:1](=O)[N:2][C:3]'
    peptide_bond_mol = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond_mol is None:
        return False, "Invalid SMARTS pattern for peptide bond"

    # Find peptide bonds
    matches = mol.GetSubstructMatches(peptide_bond_mol)
    peptide_bond_indices = []
    for match in matches:
        c_idx = match[0]
        n_idx = match[1]
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond:
            peptide_bond_indices.append(bond.GetIdx())

    # Check if there are any peptide bonds
    if not peptide_bond_indices:
        return False, "No peptide bonds found"

    # Fragment the molecule along peptide bonds
    # This will split the molecule into potential amino acid residues
    frag_mol = Chem.FragmentOnBonds(mol, peptide_bond_indices, addDummies=False)

    # Obtain the fragments
    fragments = Chem.GetMolFrags(frag_mol, asMols=True)
    num_fragments = len(fragments)

    # Check if the number of fragments (residues) is within oligopeptide range (2-20)
    if num_fragments < 2:
        return False, f"Molecule has {num_fragments} amino acid residue(s), fewer than required for an oligopeptide"
    if num_fragments > 20:
        return False, f"Molecule has {num_fragments} amino acid residues, more than allowed for an oligopeptide"

    # Define a general amino acid residue pattern
    # This pattern matches an amino group attached to an alpha carbon, which is attached to a carboxyl group
    amino_acid_pattern = Chem.MolFromSmarts('[NX3,NX4+][CX4H][CX3](=O)[O-]')
    if amino_acid_pattern is None:
        return False, "Invalid SMARTS pattern for amino acid residue"

    # Count fragments that match the amino acid pattern
    amino_acid_fragments = 0
    for frag in fragments:
        if frag.HasSubstructMatch(amino_acid_pattern):
            amino_acid_fragments += 1
        else:
            # Check for modified amino acids or terminal residues
            terminal_amino_pattern = Chem.MolFromSmarts('[NX3,NX4+][CX4H][CX3](=O)[O,N]')
            if frag.HasSubstructMatch(terminal_amino_pattern):
                amino_acid_fragments += 1

    # Validate that most fragments are amino acid residues
    if amino_acid_fragments < num_fragments:
        return False, f"Only {amino_acid_fragments} out of {num_fragments} fragments are amino acid residues"

    return True, f"Molecule has {num_fragments} amino acid residues, classified as an oligopeptide"