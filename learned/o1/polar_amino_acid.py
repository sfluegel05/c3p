"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

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

    # Identify alpha carbons connected to both amino and carboxyl groups
    alpha_carbon_pattern = Chem.MolFromSmarts('[CH1,CH2,CH3][C](N)[C](=O)O')
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon_pattern)
    if not alpha_carbons:
        return False, "No alpha carbon connected to both amino and carboxyl groups found"

    # Assume the first match is the alpha carbon
    alpha_carbon_index = alpha_carbons[0][1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_index)

    # Identify side chain atoms (exclude alpha carbon, amino group, and carboxyl group)
    side_chain_atoms = set()
    for neighbor in alpha_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != alpha_carbon_index:
            side_chain_atoms.add(neighbor.GetIdx())
            side_chain_atoms.update([atom.GetIdx() for atom in Chem.rdmolops.GetAtomNeighbors(mol, neighbor)])

    # Remove amino group and carboxyl group atoms from side chain
    amino_group = mol.GetSubstructMatch(Chem.MolFromSmarts('N'))
    carboxyl_group = mol.GetSubstructMatch(Chem.MolFromSmarts('C(=O)[O,H]'))
    exclude_atoms = set(amino_group + carboxyl_group + [alpha_carbon_index])
    side_chain_atoms = side_chain_atoms - exclude_atoms

    if not side_chain_atoms:
        return False, "No side chain found"

    # Create a side chain molecule
    side_chain = Chem.PathToSubmol(mol, side_chain_atoms)

    # Check for hydrogen bond donors and acceptors in the side chain
    hb_donor_smarts = '[N!H0,NH,H2,NH2,O!H0,OH,H2O,S!H0,SH]'
    hb_acceptor_smarts = '[O,S,N]'

    hb_donor = side_chain.GetSubstructMatch(Chem.MolFromSmarts(hb_donor_smarts))
    hb_acceptor = side_chain.GetSubstructMatch(Chem.MolFromSmarts(hb_acceptor_smarts))

    if hb_donor or hb_acceptor:
        return True, "Side chain can form hydrogen bonds"
    else:
        return False, "Side chain cannot form hydrogen bonds"