"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:27345 xanthophyll
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    Xanthophylls are oxygenated carotenoids, characterized by the presence of oxygen atoms and a conjugated double bond system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found"
    
    # Check for conjugated double bonds system indicative of carotenoids
    conjugated_systems = Chem.GetConjugatedSystems(mol)
    max_conj_length = 0
    for system in conjugated_systems:
        # Calculate the length of the conjugated system (number of atoms)
        system_length = len(system)
        if system_length > max_conj_length:
            max_conj_length = system_length
    
    # Carotenoids typically have long conjugated systems (e.g., >= 18 atoms for ~9 double bonds)
    if max_conj_length < 18:
        return False, f"Conjugated system too short ({max_conj_length} atoms)"
    
    return True, "Oxygenated carotenoid with conjugated double bond system"