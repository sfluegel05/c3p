"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: CHEBI:24400 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid must contain:
    1. A sphingoid backbone (typically with NH2 and OH groups)
    2. A fatty acid chain attached via amide bond
    3. A carbohydrate residue attached via glycosidic linkage

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible sphingoid backbone pattern
    sphingoid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4]([OX2H])[CX4H1,CX3H1]")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid backbone found"

    # Look for fatty acid chain (long carbon chain attached via amide)
    fatty_acid_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "No fatty acid chain found"

    # More comprehensive carbohydrate detection
    # Look for at least one sugar ring (pyranose or furanose)
    sugar_pattern = Chem.MolFromSmarts("[OX2][CX4H1][CX4H1][CX4H1][CX4H1][CX4H1]1[OX2][CX4H1][CX4H1][CX4H1][CX4H1]1")
    if not mol.HasSubstructMatch(sugar_pattern):
        # Try furanose pattern
        sugar_pattern = Chem.MolFromSmarts("[OX2][CX4H1][CX4H1][CX4H1][CX4H1]1[OX2][CX4H1][CX4H1][CX4H1]1")
        if not mol.HasSubstructMatch(sugar_pattern):
            return False, "No carbohydrate residue found"

    # Check that carbohydrate is attached to sphingoid via glycosidic bond
    # Look for O-C bond between sphingoid and carbohydrate
    glycosidic_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4]([OX2H])[CX4H1,CX3H1].[OX2][CX4H1][CX4H1][CX4H1][CX4H1][CX4H1]1[OX2][CX4H1][CX4H1][CX4H1][CX4H1]1")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        # Try furanose pattern
        glycosidic_pattern = Chem.MolFromSmarts("[NX3][CX4][CX4]([OX2H])[CX4H1,CX3H1].[OX2][CX4H1][CX4H1][CX4H1][CX4H1]1[OX2][CX4H1][CX4H1][CX4H1]1")
        if not mol.HasSubstructMatch(glycosidic_pattern):
            return False, "Carbohydrate not properly attached to sphingoid backbone"

    # Additional checks for typical glycosphingolipid properties
    # Molecular weight typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycosphingolipid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for glycosphingolipid"
    if o_count < 5:
        return False, "Too few oxygens for glycosphingolipid"

    return True, "Contains sphingoid backbone with fatty acid chain and carbohydrate residue"