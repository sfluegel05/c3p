"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: CHEBI:18168 galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide is a cerebroside with a galactose head group attached to a ceramide backbone.

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

    # Check for galactose head group (six-membered ring with hydroxyl groups)
    galactose_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose head group found"

    # Check for amide bond connected to long carbon chain
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Check for sphingosine/phytosphingosine backbone
    sphingosine_pattern = Chem.MolFromSmarts("[CX4][CX4]([OH])[CX4]([NH])[CX4]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine/phytosphingosine backbone found"

    # Check molecular weight - galactosylceramides are typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for galactosylceramide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for galactosylceramide"
    if o_count < 6:
        return False, "Too few oxygens for galactosylceramide"

    return True, "Contains galactose head group with ceramide backbone"