"""
Classifies: CHEBI:36500 glucosylceramide
"""
"""
Classifies: glucosylceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glucosylceramide(smiles: str):
    """
    Determines if a molecule is a glucosylceramide based on its SMILES string.
    Glucosylceramides are cerebrosides with a glucose head group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for beta-D-glucose moiety
    # Pattern for pyranose ring with correct stereochemistry
    glucose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No beta-D-glucose moiety found"

    # Check for amide linkage
    amide = Chem.MolFromSmarts("[NH][C](=O)[CH2]")
    if not mol.HasSubstructMatch(amide):
        return False, "No amide linkage found"

    # Check for sphingosine backbone with hydroxyl
    sphingosine = Chem.MolFromSmarts("[CH2]O[C]~[CH]([NH])~[CH](O)")
    if not mol.HasSubstructMatch(sphingosine):
        return False, "No sphingosine backbone found"

    # Count key atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    # Basic composition requirements
    if c_count < 20:
        return False, "Insufficient carbon atoms for glucosylceramide"
    if o_count < 7:  # Glucose has 6 oxygens + at least 1 for ceramide
        return False, "Insufficient oxygen atoms"
    if n_count != 1:
        return False, "Must have exactly one nitrogen atom"

    # Check for long alkyl chains
    alkyl_chain = "[CH2][CH2][CH2][CH2][CH2]"
    chain_matches = len(mol.GetSubstructMatches(Chem.MolFromSmarts(alkyl_chain)))
    if chain_matches < 2:
        return False, "Missing required long alkyl chains"

    # Count rotatable bonds to verify chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Insufficient rotatable bonds for required alkyl chains"

    # Verify molecular weight range (typical glucosylceramides are 600-900 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500 or mol_wt > 1200:
        return False, "Molecular weight outside typical range for glucosylceramides"

    # Check for hydroxyl groups in characteristic positions
    hydroxyl_pattern = Chem.MolFromSmarts("[CH]([OH])[CH]([NH])")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing characteristic hydroxyl groups"

    return True, "Molecule contains glucose-ceramide linkage with appropriate structural features"