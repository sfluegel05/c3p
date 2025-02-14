"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    Galactosylceramides are cerebrosides with a galactose head group and a sphingolipid backbone.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Galactose head group pattern: Galactose is typically a hexose sugar pattern that includes hydroxyl groups
    galactose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose head group found"

    # Amide bond pattern connecting a fatty acid to a sphingolipid: C(=O)N
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"
    
    # Sphingolipid backbone is characterized by the presence of an aliphatic chain and double bond(s) in certain configurations
    sphingosine_base_pattern = Chem.MolFromSmarts("[CX3](=O)N[C@@H](CO[*])[CH]=C[*]")  # adapted to capture common sphingoid bases
    if not mol.HasSubstructMatch(sphingosine_base_pattern):
        return False, "No sphingolipid backbone found"
  
    # Ensure enough oxygen for sugar and amide functionality, and enough carbon for the lipid chain
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if o_count < 6:
        return False, "Not enough oxygen atoms for galactose and amide functionality"
    if c_count < 30:
        return False, "Not enough carbon atoms for long sphingolipid and sugar chains"

    return True, "Contains galactose head group, amide linkage, and sphingolipid backbone"