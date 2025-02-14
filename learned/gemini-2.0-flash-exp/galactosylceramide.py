"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide consists of a galactose head group linked to a ceramide.
    Ceramides have a sphingosine or sphinganine backbone with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for alpha and beta galactopyranose (C1 position specified)
    alpha_galactose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    beta_galactose_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    
    #Check for a galactopyranose moiety.
    if not (mol.HasSubstructMatch(alpha_galactose_pattern) or mol.HasSubstructMatch(beta_galactose_pattern)):
            return False, "No alpha or beta galactose moiety found"

    # Relaxed SMARTS for sphingosine/phytosphingosine backbone core
    # Includes the characteristic 2-amino-1,3-diol, allowing for more variations
    sphingosine_core = Chem.MolFromSmarts("[*][C@H]([C@H](O)[C@H]([C])O)N")

    if not mol.HasSubstructMatch(sphingosine_core):
        return False, "No sphingosine/phytosphingosine backbone found"

    # Simplified glycosidic bond check: look for -C-O-C- where the two Cs are in galactose and sphingosine respectively
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C]([OH1])([C@H]([C@H]([C@H]([C@H](O[*])CO)O)O)O)[O][C]([C@H](O)[C@H]([C])O)[N]")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No glycosidic bond present between galactose and ceramide"
        
    # Improved fatty acid chain check (include C=O)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX3](=[OX1])")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"

    return True, "Contains a galactose moiety linked to a ceramide backbone."