"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid containing one or more sialic acids 
    linked to an oligosaccharide chain attached to a ceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight - gangliosides are large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for ganglioside"

    # Look for ceramide backbone components with more flexible patterns
    # Multiple possible sphingosine patterns to catch variations
    sphingosine_patterns = [
        Chem.MolFromSmarts("[CH2,CH][CH2,CH][CH]([OH,O])[CH]([NH,N])[CH2]O"), # Basic pattern
        Chem.MolFromSmarts("C=CC[CH]([OH,O])[CH]([NH,N])[CH2]O"),  # With double bond
        Chem.MolFromSmarts("[CH2,CH][CH2,CH][CH]([OH,O])[CH]([NH,N])CO")  # Alternative connection
    ]
    
    # Multiple possible acyl patterns
    acyl_patterns = [
        Chem.MolFromSmarts("[NH,N][C](=O)[CH2][CH2,CH3]"),  # Basic pattern
        Chem.MolFromSmarts("[NH,N]C(=O)CC"),  # Shorter version
        Chem.MolFromSmarts("[NH,N]C(=O)CCCCC")  # Longer explicit version
    ]
    
    # Check if any combination of patterns matches
    has_sphingosine = any(mol.HasSubstructMatch(pattern) for pattern in sphingosine_patterns)
    has_acyl = any(mol.HasSubstructMatch(pattern) for pattern in acyl_patterns)
    
    if not (has_sphingosine and has_acyl):
        return False, "No ceramide backbone found"

    # Look for sialic acid (Neu5Ac) patterns - multiple possibilities
    sialic_patterns = [
        Chem.MolFromSmarts("[C,c](=O)[O;H,-][C,c]1[C,c][C,c]([NH,N][C](=O))[C,c][C,c][C,c]1"),  # Basic pattern
        Chem.MolFromSmarts("C(=O)[O;H,-]C1CC(NC=O)CC(O)C1"),  # More specific pattern
        Chem.MolFromSmarts("C(=O)[O]C1C(O)C(N)CC(O)C1")  # Alternative pattern
    ]
    
    has_sialic = False
    for pattern in sialic_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            has_sialic = True
            break
            
    if not has_sialic:
        return False, "No sialic acid (Neu5Ac) found"

    # Look for oligosaccharide patterns
    sugar_patterns = [
        Chem.MolFromSmarts("[C,c]1[O,o][C,c][C,c][C,c][C,c]1"),  # Basic pyranose
        Chem.MolFromSmarts("C1OCCC(O)C1"),  # Simplified sugar
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1")  # More specific pattern
    ]
    
    sugar_count = 0
    for pattern in sugar_patterns:
        if pattern is not None:
            sugar_count += len(mol.GetSubstructMatches(pattern))
    
    if sugar_count < 3:
        return False, "Insufficient oligosaccharide chain"

    # Check for glycosidic linkages
    glycosidic_pattern = Chem.MolFromSmarts("[C,c][O,o][C,c]")
    if glycosidic_pattern is not None:
        glycosidic_count = len(mol.GetSubstructMatches(glycosidic_pattern))
        if glycosidic_count < 4:
            return False, "Insufficient glycosidic linkages"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 30:
        return False, "Too few carbons for ganglioside structure"
    if o_count < 15:
        return False, "Too few oxygens for ganglioside structure"
    if n_count < 2:
        return False, "Too few nitrogens for ganglioside structure"

    return True, "Contains ceramide backbone, sialic acid(s), and oligosaccharide chain with correct linkages"