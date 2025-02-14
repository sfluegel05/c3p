"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid is an amino acid whose side chain is capable of forming one or more hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the amino acid backbone pattern
    backbone_pattern = Chem.MolFromSmarts("[NX3,NH2,NH3+][CX4H]([*])[CX3](=O)[O-,OH]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Molecule does not have a standard amino acid backbone"

    # Find the match for the backbone to identify backbone atoms
    backbone_match = mol.GetSubstructMatch(backbone_pattern)
    backbone_atoms = set(backbone_match)

    # Identify side chain atoms (atoms not in backbone)
    side_chain_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in backbone_atoms:
            side_chain_atoms.add(atom.GetIdx())

    if not side_chain_atoms:
        return False, "No side chain found"

    # Check for hydrogen bond donors and acceptors in side chain
    # Functional groups capable of hydrogen bonding: -OH, -SH, -NH2, -CONH2, imidazole, phenol, carboxylic acids, etc.
    hbond_donor_acceptor_smarts = [
        "[OX2H]",          # Hydroxyl group
        "[SX2H]",          # Thiol group
        "[NX3H2]",         # Primary amine
        "[NX3H][CX3](=O)[#6]",  # Amide group
        "c1nccc1",         # Imidazole ring
        "c1cc([OX2H])ccc1",# Phenol group
        "[CX3](=O)[OX1H0-]"# Carboxylate group
    ]
    side_chain_mol = Chem.PathToSubmol(mol, list(side_chain_atoms))
    for smarts in hbond_donor_acceptor_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if side_chain_mol.HasSubstructMatch(pattern):
            return True, "Side chain can form hydrogen bonds"

    return False, "Side chain cannot form hydrogen bonds"