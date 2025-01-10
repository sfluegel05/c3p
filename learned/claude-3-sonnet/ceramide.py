"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: ceramides (N-acyl-sphingoid bases)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide bond (-NH-C(=O)-)
    amide_pattern = Chem.MolFromSmarts("[NX3H]-[CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"
    
    # Count amide bonds - ceramides should have exactly one
    amide_matches = len(mol.GetSubstructMatches(amide_pattern))
    if amide_matches != 1:
        return False, f"Found {amide_matches} amide bonds, need exactly 1"

    # Look for hydroxyl groups
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 2:
        return False, f"Found only {oh_matches} hydroxyl groups, need at least 2"

    # Check for long carbon chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long carbon chain found"

    # Count carbons and check molecular weight
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Only {c_count} carbons found, need at least 20"

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for typical ceramide"

    # Look for primary alcohol (-CH2-OH) near amide
    primary_alcohol_pattern = Chem.MolFromSmarts("[NX3]-[CX4]-[CX4]-[OX2H]")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol found near amide group"

    # Additional checks for common modifications
    has_sugar = False
    has_phosphate = False
    has_sulfate = False
    
    # Check for glycosylation (presence of cyclic acetal/ketal)
    sugar_pattern = Chem.MolFromSmarts("[OX2]1[CX4][CX4][CX4][CX4][CX4]1")
    if mol.HasSubstructMatch(sugar_pattern):
        has_sugar = True
        
    # Check for phosphate
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2H,OX1-])[OX2H,OX1-]")
    if mol.HasSubstructMatch(phosphate_pattern):
        has_phosphate = True
        
    # Check for sulfate
    sulfate_pattern = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX2H,OX1-]")
    if mol.HasSubstructMatch(sulfate_pattern):
        has_sulfate = True

    # Build classification reason
    mods = []
    if has_sugar:
        mods.append("glycosylated")
    if has_phosphate:
        mods.append("phosphorylated")
    if has_sulfate:
        mods.append("sulfated")
    
    mod_str = " and ".join(mods)
    if mod_str:
        reason = f"Ceramide ({mod_str})"
    else:
        reason = "Basic ceramide structure"

    return True, f"{reason} with long-chain base, amide-linked fatty acid, and characteristic hydroxylation pattern"