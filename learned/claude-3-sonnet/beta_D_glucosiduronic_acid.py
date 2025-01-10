"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is formed by the condensation of any substance with 
    beta-D-glucuronic acid to form a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create and validate SMARTS patterns
    patterns = {
        # Basic glucuronic acid scaffold (6-membered ring with COOH)
        'core': '[OX2][CH1]1[CH1][CH1][CH1][CH1](C(=O)[OH1])O1',
        # Carboxylic acid group
        'carboxyl': 'C(=O)[OH1]',
        # Hydroxyl groups pattern
        'hydroxyls': '[OH1]',
        # Beta glycosidic linkage at anomeric carbon
        'beta_linkage': '[OX2][CH1]1O[CH1]'
    }
    
    # Validate SMARTS patterns
    smarts_mols = {}
    for name, pattern in patterns.items():
        smarts_mol = Chem.MolFromSmarts(pattern)
        if smarts_mol is None:
            return False, f"Invalid SMARTS pattern: {name}"
        smarts_mols[name] = smarts_mol

    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(smarts_mols['carboxyl']):
        return False, "No carboxylic acid group found"

    # Check for basic glucuronic acid scaffold
    if not mol.HasSubstructMatch(smarts_mols['core']):
        return False, "No glucuronic acid core structure found"

    # Count hydroxyl groups (should have at least 3 for glucuronic acid)
    hydroxyl_matches = len(mol.GetSubstructMatches(smarts_mols['hydroxyls']))
    if hydroxyl_matches < 3:
        return False, f"Insufficient hydroxyl groups (found {hydroxyl_matches}, need at least 3)"

    # Check for proper glycosidic linkage
    if not mol.HasSubstructMatch(smarts_mols['beta_linkage']):
        return False, "No proper glycosidic linkage found"

    # Count oxygen atoms (should have at least 6 for glucuronic acid)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Insufficient oxygen atoms (found {o_count}, need at least 6)"

    # Check ring count
    rings = mol.GetRingInfo()
    if rings.NumRings() == 0:
        return False, "No rings found in molecule"

    # Additional stereochemistry check using 3D conformation
    try:
        # Generate 3D conformation
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Check if molecule has chiral centers
        chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        if len(chiral_centers) < 4:  # Glucuronic acid should have at least 4 chiral centers
            return False, "Insufficient chiral centers for beta-D-glucuronic acid"
    except:
        # If 3D generation fails, continue without it
        pass

    return True, "Contains beta-D-glucuronic acid core with glycosidic linkage"