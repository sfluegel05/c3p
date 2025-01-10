"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
from rdkit import Chem

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Dihydroagarofuran sesquiterpenoids are sesquiterpenoids with a specific bicyclic skeleton.

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

    # Updated SMARTS pattern for dihydroagarofuran skeleton (hypothetical example)
    dihydroagarofuran_pattern = Chem.MolFromSmarts("[C@]12(CO[C@H]3OC(=O)[C@@H]([C@H]1[C@H]4OC[C@H]2O4)[C@@H]3C(C)=O)") 
    
    if dihydroagarofuran_pattern is None or not mol.HasSubstructMatch(dihydroagarofuran_pattern):
        return False, "No dihydroagarofuran skeleton found"

    # Check for sesquiterpenoid characteristics
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15:
        return False, f"Too few carbons for a sesquiterpenoid, found {c_count} carbon atoms"

    # Check for presence of typical functional groups (e.g., ester or ketone groups)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1,C]")  # esters and acids
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester or acid groups found, which are common in this class"

    return True, "Contains dihydroagarofuran skeleton and sesquiterpenoid characteristics"

# Example testing
# Expected result based on known classification of the input
print(is_dihydroagarofuran_sesquiterpenoid("C1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(C)=O)[C@@H](O)[C@H](C[C@@H](C)[C@@]3(COC1(C)C)C)[C@@H]2OC(=O)c1ccccc1)OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O"))