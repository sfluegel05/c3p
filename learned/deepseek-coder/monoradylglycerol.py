"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: CHEBI:75568 monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol backbone with a single acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 hydroxyls and 1 substituent)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OH])[CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester (acyl), alkyl, or alk-1-enyl substituent
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")  # Acyl (ester)
    alkyl_pattern = Chem.MolFromSmarts("[CX4]")  # Alkyl
    alkenyl_pattern = Chem.MolFromSmarts("[CX3]=[CX3]")  # Alk-1-enyl

    # Count the number of substituents
    substituent_count = 0
    if mol.HasSubstructMatch(ester_pattern):
        substituent_count += 1
    if mol.HasSubstructMatch(alkyl_pattern):
        substituent_count += 1
    if mol.HasSubstructMatch(alkenyl_pattern):
        substituent_count += 1

    # Check if there is exactly one substituent
    if substituent_count != 1:
        return False, f"Found {substituent_count} substituents, need exactly 1"

    # Check molecular weight - monoradylglycerols typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for monoradylglycerol"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 6:
        return False, "Too few carbons for monoradylglycerol"
    if o_count < 3:
        return False, "Must have at least 3 oxygens (glycerol backbone)"

    return True, "Contains glycerol backbone with a single acyl, alkyl, or alk-1-enyl substituent"