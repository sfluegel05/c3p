"""
Classifies: CHEBI:60027 polymer
"""
from rdkit import Chem
from rdkit.Chem import rdmolops, rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    Focuses on identifying extensive chains, repeating units, and significant molecular size—common polymer features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is identified as a polymer-like structure, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string to create a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define polymer-like repeating unit patterns with more representative recurrence
    patterns = [
        Chem.MolFromSmarts("[C](-[C,C])(-[C,C])[C]"),  # generic polymer backbone
        Chem.MolFromSmarts("C-C-C-C"),  # flexible aliphatic chain
        Chem.MolFromSmarts("c1ccccc1"),  # aromatic polymer (e.g., polystyrenes)
        Chem.MolFromSmarts("[N,O]-[C]-[O,N]"),  # amide/ester linkages common in polyamides, polyesters
    ]

    # Check matches for polymer-like chains
    total_matches = sum(len(mol.GetSubstructMatches(pattern)) for pattern in patterns)

    # Refined logic: multiple criteria check
    num_atoms = mol.GetNumAtoms()
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    if total_matches > 10:
        return True, f"Contains significant repeating substructural patterns: found {total_matches}"
    
    if num_atoms > 100 and total_matches > 5:
        return True, f"Large molecular size with polymer-linked repeating units, {num_atoms} atoms"
    
    if rotatable_bonds > 50:
        return True, f"High flexibility and length hinting at polymeric property, {rotatable_bonds} rotatable bonds"
    
    return False, "No distinct polymer characteristics detected"