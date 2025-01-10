"""
Classifies: CHEBI:16219 cucurbitacin
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane,
    typically featuring extensive hydroxylation and distinctive carbonyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least 4 rings (approximation for tetracyclic structure)
    if len(Chem.GetSymmSSSR(mol)) < 4:
        return False, "Less than 4 rings, does not meet tetracyclic criteria"
    
    # Check for hydroxyl groups using SMARTS [#6][OH]
    hydroxyl_smarts = Chem.MolFromSmarts('[CX4][OH]')
    if not mol.HasSubstructMatch(hydroxyl_smarts):
        return False, "Missing expected hydroxyl groups"
    
    # Check for carbonyl groups using SMARTS C=O
    carbonyl_smarts = Chem.MolFromSmarts('C=O')
    if not mol.HasSubstructMatch(carbonyl_smarts):
        return False, "Missing expected carbonyl groups"

    # Count total carbon atoms, cucurbitacins typically have 30 carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, "Too few carbons to represent typical cucurbitacins"
    
    # Provide a more detailed structural match, targeting the cucurbitane backbone
    cucurbitane_smarts = Chem.MolFromSmarts('C1CCC2C3C1C(=O)C=C4C3C(C2)[C@@]5(CC[C@H](O)[C@@H](O)C5=O)C4')
    if not mol.HasSubstructMatch(cucurbitane_smarts):
        return False, "Does not contain the cucurbitane backbone"
    
    return True, "This molecule matches cucurbitacin structural features typical of tetracyclic triterpenoids"