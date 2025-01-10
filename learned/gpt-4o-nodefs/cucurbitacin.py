"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are highly oxygenated tetracyclic triterpenes with specific structural features.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Defining a basic pattern for cucurbitacins
    tetracyclic_pattern = Chem.MolFromSmarts('[C@]12CC[C@@H]3[C@@]1(C)C(=O)C[C@]4(C)[C@@H]3CC[C@@]4([H])[C@]2(C)[C@@](C)(O)[C]=[O]')
    if not mol.HasSubstructMatch(tetracyclic_pattern):
        return False, "Does not match tetracyclic skeleton pattern"

    # Additional features like high oxygenation
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Check if there's a reasonable ratio of oxygen to carbons typical of cucurbitacins
    if oxygen_count / carbon_count < 0.2:  # This is arbitrary; adjust based on characteristic ratios in cucurbitacins
        return False, "Insufficient oxygenation characteristic of cucurbitacins"

    # Further features and checks may be added for higher specificity
    return True, "Contains features characteristic of cucurbitacins"