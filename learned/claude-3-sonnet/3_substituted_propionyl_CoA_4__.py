"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    Must have propionyl base structure substituted at 3-position and exactly 4 negative charges.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 4 negative charges
    charge_pattern = Chem.MolFromSmarts("[O-]")
    charge_matches = mol.GetSubstructMatches(charge_pattern)
    if len(charge_matches) != 4:
        return False, f"Found {len(charge_matches)} negative charges, need exactly 4"

    # Check for complete CoA core structure
    coa_core = Chem.MolFromSmarts("[$(N1C=NC=2C(N)=NC=NC2=1)]" # Adenine
                                 "[$(OC1C(O)C(COP([O-])([O-])=O)OC1)]" # Ribose + phosphate
                                 "[$(OP(=O)([O-])OP(=O)([O-])OCC(C)(C)[CH](O))]" # Diphosphate + pantothenic acid
                                 "[$(C(=O)NCCC(=O)NCCSC(=O))]") # Pantetheine arm + thioester
    if not mol.HasSubstructMatch(coa_core):
        return False, "Missing or incomplete CoA core structure"

    # Check for 3-substituted propionyl structure (C-C-C with substitution at C3)
    # The pattern looks for: S-C(=O)-CH2-CH(R)-R where R is not H
    propionyl_pattern = Chem.MolFromSmarts("SC(=O)CC([!H])([!H])")
    if not mol.HasSubstructMatch(propionyl_pattern):
        return False, "No 3-substituted propionyl structure found"

    # Exclude specific patterns that would indicate different CoA derivatives
    exclude_patterns = [
        ("SC(=O)C(O)", "2-hydroxy"), # 2-hydroxy CoAs
        ("SC(=O)CC(=O)", "3-oxo"), # 3-oxo CoAs
        ("SC(=O)CCO", "hydroxy"), # hydroxy CoAs
        ("SC(=O)C(=O)", "2-oxo"), # 2-oxo CoAs
    ]
    
    for pattern, name in exclude_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Molecule is a {name}-CoA derivative"

    # Count key atoms to verify overall composition
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    
    if p_count != 3:
        return False, f"Must have exactly 3 phosphorus atoms, found {p_count}"
    if s_count != 1:
        return False, f"Must have exactly 1 sulfur atom, found {s_count}"

    return True, "Contains CoA(4-) moiety with 3-substituted propionyl structure"