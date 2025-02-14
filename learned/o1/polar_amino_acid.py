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

    # Define the amino acid backbone pattern (inclusive of zwitterions)
    backbone_pattern = Chem.MolFromSmarts("[N;H1,H2;+0;+1][C@@,H](C)[C](=O)[O;H1,-1]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Molecule does not have a standard amino acid backbone"

    # Find backbone atoms
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if not backbone_matches:
        return False, "Amino acid backbone not found"
    backbone_atoms = set(backbone_matches[0])

    # Identify side chain atoms (atoms not in backbone)
    side_chain_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() not in backbone_atoms]
    if not side_chain_atoms:
        return False, "No side chain found"

    # Define SMARTS patterns for hydrogen bonding groups
    hbond_donor_acceptor_smarts = [
        "[OX1H0-,OX2H1]",        # Hydroxyl group and oxyanions
        "[SX2H]",                # Thiol group
        "[NX3H2,NX3H1,NX4+]",    # Primary and secondary amines, protonated amines
        "[NX3][CX3](=O)[NX3]",   # Amide group
        "c1ncn[cH]1",            # Imidazole group in histidine
        "c1ccc(O)cc1",           # Phenol group in tyrosine
        "C(=O)[O-,OH]"           # Carboxylate and carboxylic acid
    ]

    # Check for hydrogen bonding groups in the side chain
    for atom in side_chain_atoms:
        for smarts in hbond_donor_acceptor_smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if mol.HasSubstructMatch(pattern, atomIndices=[atom.GetIdx()]):
                return True, "Side chain contains hydrogen bonding groups"

    return False, "Side chain does not contain hydrogen bonding groups"