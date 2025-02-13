"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: galactosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide has a galactose sugar connected to a ceramide 
    (sphingosine + fatty acid) via a glycosidic bond.

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

    # Look for pyranose ring (6-membered sugar ring with oxygen)
    # This pattern matches the basic ring structure of galactose
    sugar_ring = Chem.MolFromSmarts("[C]-1-[C]-[C]-[C]-[C]-[O]-1")
    if not mol.HasSubstructMatch(sugar_ring):
        return False, "No sugar ring found"

    # Look for characteristic galactose hydroxyl/substitution pattern
    # Allow for both OH and OS (sulfated) groups
    galactose_pattern = Chem.MolFromSmarts("[C]-1-[C]([O])-[C]([O])-[C]([O])-[C]([CH2][O])-[O]-1")
    if not mol.HasSubstructMatch(galactose_pattern):
        return False, "No galactose sugar pattern found"

    # Look for amide bond characteristic of ceramide
    amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")
    if not mol.HasSubstructMatches(amide_pattern):
        return False, "No amide bond found"

    # Look for sphingosine/ceramide backbone with characteristic OH groups
    sphingosine_pattern = Chem.MolFromSmarts("[CH2][O]-[CH]-[CH]([O])-[CH]([NH])-[CH2]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine/ceramide backbone found"

    # Check for long carbon chains (at least 12 carbons in a row)
    long_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "Missing required long hydrocarbon chain"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for galactosylceramide"
    if o_count < 6:
        return False, "Too few oxygens for galactosylceramide"
    if n_count != 1:
        return False, "Must have exactly one nitrogen (amide bond)"

    # Check for optional sulfate group
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O,OH]")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    base_reason = "Contains galactose connected to ceramide via glycosidic bond"
    if has_sulfate:
        return True, base_reason + " with sulfate modification"
    return True, base_reason