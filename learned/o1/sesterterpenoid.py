"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from sesterterpenes composed of five isoprene units (C25 skeleton),
    possibly modified by rearrangement or removal of small groups (e.g., methyl groups).
    They typically contain ring structures and methyl branching.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 70:
        return False, f"Carbon count is {c_count}, not typical for sesterterpenoids"

    # Check for ring structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structures found; molecule is likely acyclic"

    # Calculate the length of the longest aliphatic chain
    longest_chain = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            paths = Chem.rdmolops.FindAllPathsOfLengthN(mol, 1, useBonds=True)
            for n in range(1, c_count + 1):
                paths = Chem.rdmolops.FindAllPathsOfLengthN(mol, n, useBonds=True)
                for path in paths:
                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path):
                        chain_length = len(path)
                        if chain_length > longest_chain:
                            longest_chain = chain_length
    if longest_chain > 15:
        return False, f"Longest aliphatic chain is {longest_chain}, indicates a linear molecule"

    # Check for methyl branching (tertiary carbons)
    tertiary_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and len(atom.GetNeighbors()) == 3]
    if len(tertiary_carbons) < 3:
        return False, f"Found {len(tertiary_carbons)} tertiary carbons; terpenoids typically have methyl branching"
    
    # Check for oxygen-containing functional groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; terpenoids are often oxygenated"

    # Check for typical terpenoid functional groups
    # e.g., hydroxyl, carbonyl, ester, ether
    functional_groups = [
        Chem.MolFromSmarts('[OX2H]'),     # Hydroxyl group
        Chem.MolFromSmarts('C=O'),        # Carbonyl group
        Chem.MolFromSmarts('COC'),        # Ether group
        Chem.MolFromSmarts('OC=O'),       # Ester group
    ]
    fg_found = False
    for fg in functional_groups:
        if mol.HasSubstructMatch(fg):
            fg_found = True
            break
    if not fg_found:
        return False, "No typical terpenoid functional groups found"

    return True, "Molecule meets the refined criteria for a sesterterpenoid"