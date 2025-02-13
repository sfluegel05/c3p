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

    # Look for galactose pattern - 6-membered ring with specific OH pattern
    # More specific pattern for galactose
    galactose_pattern = Chem.MolFromSmarts("[CH2][OH]-[CH]1-[CH]([OH])-[CH]([OH])-[CH]([OH])-[CH]([OH])-O1")
    if not mol.HasSubstructMatch(galactose_pattern):
        # Try alternative pattern that might match sulfated versions
        alt_galactose = Chem.MolFromSmarts("[CH2][OH]-[CH]1-[CH]([OH])-[CH]([O])-[CH]([OH])-[CH]-O1")
        if not mol.HasSubstructMatch(alt_galactose):
            return False, "No galactose sugar found"

    # Look for amide bond (-NH-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    # Look for sphingosine/ceramide backbone
    # Pattern matches the characteristic OH-NH-OH region with attached chains
    sphingosine_pattern = Chem.MolFromSmarts("[CH2]O[CH]-[CH]([OH])-[CH]([NH])-[CH2]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine/ceramide backbone found"

    # Check for long carbon chains
    long_chain = Chem.MolFromSmarts("[CH2][CH2][CH2][CH2][CH2][CH2][CH2]")
    chain_matches = len(mol.GetSubstructMatches(long_chain))
    if chain_matches < 2:
        return False, "Missing required long hydrocarbon chains"

    # Count atoms to verify molecule size
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
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[OH]")
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    base_reason = "Contains galactose connected to ceramide via glycosidic bond"
    if has_sulfate:
        return True, base_reason + " with sulfate modification"
    return True, base_reason