"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is characterized by having an ester linkage with an 
    8-carbon chain forming part of the ester functional group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General ester pattern: a carbonyl (C=O) with a single-bonded oxygen attached
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Identify ester groups and validate the presence of an octanoate component
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for match in ester_matches:
        # Validate if a connected chain is 8 carbons long (to form octanoate)
        atom_idx = match[2]  # Atom where the ester oxygen is attached
        octanoate_frag = Chem.rdmolops.FragmentOnBonds(mol, [atom_idx], addDummies=True)
        for frag in Chem.GetMolFrags(octanoate_frag):
            fragment_mol = Chem.MolFragmentToSmiles(octanoate_frag, atomsToUse=frag)
            fragment_mol = Chem.MolFromSmiles(fragment_mol)
            carbon_chain_length = sum(1 for atom in fragment_mol.GetAtoms() if atom.GetAtomicNum() == 6)
            if carbon_chain_length == 8:
                return True, "Contains octanoate ester group"
    
    return False, "No octanoate ester group found"