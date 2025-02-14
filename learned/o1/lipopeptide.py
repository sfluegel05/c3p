"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: lipopeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with an attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify peptide bonds (amide bonds between amino acids)
    peptide_bond_smarts = "[C;!R](=O)N"  # Non-ring carbonyl carbon bonded to nitrogen
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)
    if num_peptide_bonds < 2:
        return False, f"Insufficient peptide bonds found ({num_peptide_bonds} found, need at least 2)"

    # Identify amino acid residues (backbone pattern)
    amino_acid_smarts = "[N;$([N;H1,H2][C;!R]);!R][C;!R](=O)"  # N-C(=O)
    amino_acid_pattern = Chem.MolFromSmarts(amino_acid_smarts)
    amino_acids = mol.GetSubstructMatches(amino_acid_pattern)
    num_amino_acids = len(amino_acids)
    if num_amino_acids < 2:
        return False, f"Insufficient amino acid residues found ({num_amino_acids} found, need at least 2)"

    # Identify long aliphatic chains (lipid chains)
    # We will look for chains of 8 or more carbons not in rings and not attached to heteroatoms
    max_chain_length = 0
    for chain in Chem.FindAllPathsOfLengthN(mol, 8, useBonds=True):
        length = 0
        for idx in chain:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.IsInRing():
                break
            hetero_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 6]
            if hetero_neighbors:
                break
            length += 1
        if length == 8:
            max_chain_length = length
            break  # Found a chain of sufficient length

    if max_chain_length < 8:
        return False, f"No lipid chain found (longest aliphatic chain is {max_chain_length} carbons, need at least 8)"

    # Check connectivity between lipid chain and peptide
    # Identify the lipid chain atoms
    lipid_smarts = "[C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R]"  # Chain of 8 aliphatic carbons
    lipid_pattern = Chem.MolFromSmarts(lipid_smarts)
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)

    if not lipid_matches:
        return False, "No lipid chain matching pattern found"

    # Check if any lipid chain is connected to the peptide
    lipid_atoms = set()
    for match in lipid_matches:
        lipid_atoms.update(match)

    peptide_atoms = set()
    for bond in peptide_bonds:
        peptide_atoms.update(bond)

    # Find if there's a path between any lipid atom and any peptide atom
    is_connected = False
    for lipid_atom in lipid_atoms:
        for peptide_atom in peptide_atoms:
            if Chem.rdmolops.HasPath(mol, lipid_atom, peptide_atom):
                is_connected = True
                break
        if is_connected:
            break

    if not is_connected:
        return False, "Lipid chain is not connected to peptide moiety"

    # Additional filter: check for ester or amide linkage between lipid and peptide
    linkage_found = False
    ester_or_amide_smarts = "[C;!R](=O)[O,N][C;!R]"  # Ester or amide linkage
    linkage_pattern = Chem.MolFromSmarts(ester_or_amide_smarts)
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)
    for match in linkage_matches:
        if match[0] in lipid_atoms and match[-1] in peptide_atoms:
            linkage_found = True
            break
        if match[-1] in lipid_atoms and match[0] in peptide_atoms:
            linkage_found = True
            break

    if not linkage_found:
        return False, "No linkage found between lipid chain and peptide moiety"

    return True, "Contains both peptide bonds and lipid chain connected by appropriate linkage indicative of a lipopeptide"