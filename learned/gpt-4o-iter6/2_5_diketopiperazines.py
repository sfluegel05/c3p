"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a 2,5-diketopiperazine with allowance for stereochemistry and substituents
    diketopiperazine_core = Chem.MolFromSmarts("N1C(=O)[C@H]2NC(=O)[C@@H]2C1")  # Incorporates stereochemistry and ring closure
    
    if mol.HasSubstructMatch(diketopiperazine_core):
        # Ensuring the matched substructure is indeed piperazine-2,5-dione and not a similar motif
        # Ensure piperazine ring by counting exactly 6 atoms in the core structure
        if sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic()) == 0:  # Piperazines are not aromatic
            core_atoms = mol.GetSubstructMatch(diketopiperazine_core)
            if core_atoms and len(core_atoms) == 6:
                return True, "Contains a piperazine-2,5-dione skeleton"
            else:
                return False, "Molecule fails 6-membered core match criterion"

    return False, "Does not contain a piperazine-2,5-dione skeleton"