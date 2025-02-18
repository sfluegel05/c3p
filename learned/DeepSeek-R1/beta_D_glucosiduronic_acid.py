"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
"""
Classifies: CHEBI:16199 beta-D-glucosiduronic acid
"""
from rdkit import Chem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    Must contain a beta-D-glucuronic acid moiety linked via glycosidic bond.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Enhanced SMARTS pattern for beta-D-glucuronic acid with glycosidic bond:
    # 1. [O;X2] = glycosidic oxygen (must have two bonds)
    # 2. Correct stereochemistry at all chiral centers (beta-D configuration)
    # 3. Carboxylic acid group at C5 (C(=O)O)
    glucuronic_pattern = Chem.MolFromSmarts(
        "[O;X2][C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H](C(=O)O)O1"
    )
    
    # Verify pattern match and check for at least one glycosidic bond
    if mol.HasSubstructMatch(glucuronic_pattern):
        # Additional check for actual glycosidic linkage (O connected to non-H)
        for match in mol.GetSubstructMatches(glucuronic_pattern):
            glycosidic_oxygen_idx = match[0]
            atom = mol.GetAtomWithIdx(glycosidic_oxygen_idx)
            # Ensure oxygen is connected to another atom (not just part of the ring)
            if atom.GetDegree() >= 2:
                return True, "Contains beta-D-glucuronic acid with glycosidic bond"
        
    return False, "No beta-D-glucuronic acid glycoside detected"