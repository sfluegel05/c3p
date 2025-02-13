"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid with a carbohydrate attached via glycosidic linkage 
    to a sphingoid/ceramide base.

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

    # Look for ceramide/sphingoid base patterns
    # NH-C(=O) amide group
    amide_pattern = Chem.MolFromSmarts("[NH][C;X3](=[O;X1])")
    # Long carbon chain after amide
    chain_pattern = Chem.MolFromSmarts("[C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]")
    
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found (required for ceramide/sphingoid base)"
    
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long carbon chain found (required for ceramide/sphingoid base)"

    # Look for sugar patterns
    # Pyranose ring with multiple OH groups
    sugar_pattern = Chem.MolFromSmarts("[C;R1][O;R1][C;R1][C;R1][C;R1][C;R1]")
    # Multiple OH groups
    oh_pattern = Chem.MolFromSmarts("[OH]")
    
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar ring found"
    
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:
        return False, "Not enough hydroxyl groups for sugar moiety"

    # Look for glycosidic linkage (C-O-C)
    glycosidic_pattern = Chem.MolFromSmarts("[C;R1][O;X2][C;X4]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Additional check for sphingoid characteristics
    # Look for OH group near amide and long chain
    sphingoid_oh_pattern = Chem.MolFromSmarts("[NH]C([C;X4])[C;X4][OH]")
    if not mol.HasSubstructMatch(sphingoid_oh_pattern):
        return False, "Missing characteristic sphingoid OH group"

    # Count carbons and oxygens to ensure reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for glycosphingolipid"
    if o_count < 6:
        return False, "Too few oxygens for glycosphingolipid"

    return True, "Contains carbohydrate residue attached via glycosidic linkage to sphingoid/ceramide base"