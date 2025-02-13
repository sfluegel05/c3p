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

    # Look for ceramide backbone components
    # More flexible sphingosine pattern - long chain base with OH and NH
    sphingosine_pattern = Chem.MolFromSmarts("[CH2,CH]-[CH2,CH]-[CH](-[OH])-[CH](-[NH])-[CH2]-O")
    # Acyl chain attached to NH
    acyl_pattern = Chem.MolFromSmarts("[NH]-[C](=O)-[CH2]-[CH2]")
    
    if not (mol.HasSubstructMatch(sphingosine_pattern) and mol.HasSubstructMatch(acyl_pattern)):
        return False, "No ceramide backbone found"

    # Look for sialic acid (Neu5Ac) pattern - more flexible version
    # Core structure with carboxyl and N-acetyl
    sialic_pattern = Chem.MolFromSmarts("[C,c](=O)[O;H,-]-[C,c]1-[C,c]-[C,c](-[NH]-[C](=O))-[C,c]-[C,c]-[C,c]1")
    sialic_matches = mol.GetSubstructMatches(sialic_pattern)
    if len(sialic_matches) < 1:
        return False, "No sialic acid (Neu5Ac) found"

    # Look for oligosaccharide pattern - more flexible sugar detection
    sugar_pattern = Chem.MolFromSmarts("[C,c]1[O,o][C,c][C,c][C,c][C,c]1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 3:
        return False, "Insufficient oligosaccharide chain"

    # Check for glycosidic linkages between sugars
    glycosidic_pattern = Chem.MolFromSmarts("[C,c]-[O,o]-[C,c]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 4:
        return False, "Insufficient glycosidic linkages"

    # Count key atoms to verify overall composition
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 30:
        return False, "Too few carbons for ganglioside structure"
    if o_count < 15:
        return False, "Too few oxygens for ganglioside structure"
    if n_count < 2:
        return False, "Too few nitrogens for ganglioside structure"

    # Verify presence of long chain fatty acid - more flexible pattern
    fatty_acid_pattern = Chem.MolFromSmarts("[CH2,CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "Missing long chain fatty acid"

    return True, "Contains ceramide backbone, sialic acid(s), and oligosaccharide chain with correct linkages"