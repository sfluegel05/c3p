"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Dihydroagarofuran sesquiterpenoids are sesquiterpenoids with a dihydroagarofuran skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for dihydroagarofuran skeleton
    # Replace with the actual SMARTS pattern for dihydroagarofuran
    dihydroagarofuran_pattern = Chem.MolFromSmarts("C1[C@]23COC(=O)C[C@H]2C4[C@@H]1CO[C@H]34")  # Example pattern
    
    if dihydroagarofuran_pattern is None or not mol.HasSubstructMatch(dihydroagarofuran_pattern):
        return False, "No dihydroagarofuran skeleton found"

    # Check for sesquiterpenoid characteristics
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbons for a sesquiterpenoid, found {c_count} carbon atoms"

    # The exact check might vary depending on specifics, add additional characteristics if known.

    return True, "Contains dihydroagarofuran skeleton and has sesquiterpenoid characteristics"

# Example testing
# Expected result based on known classification of the input
print(is_dihydroagarofuran_sesquiterpenoid("C1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(C)=O)[C@@H](O)[C@H](C[C@@H](C)[C@@]3(COC1(C)C)C)[C@@H]2OC(=O)c1ccccc1)OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O"))