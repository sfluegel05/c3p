"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
"""
Classifies: CHEBI:35943 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is a ceramide consisting of sphinganine with one of the amino hydrogens
    substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphinganine backbone pattern (C-C-N-C-C-O)
    sphinganine_pattern = Chem.MolFromSmarts("[NH1X3,NH2X2,NH1X3+]-[CH2X4]-[CH1X4]-[CH2X4]-[CH2X4]-[OX2H1]")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"

    # Look for fatty acyl group attached to N
    fatty_acyl_pattern = Chem.MolFromSmarts("[NH1X3,NH2X2,NH1X3+]-[CX3](=[OX1])-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)
    if not fatty_acyl_matches:
        return False, "No fatty acyl group attached to N"

    # Check for long carbon chains (fatty acids typically > 10 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Fatty acid chain too short"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - N-acylsphinganines typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for N-acylsphinganine"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for N-acylsphinganine"
    if o_count < 4:
        return False, "Too few oxygens for N-acylsphinganine"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen"

    return True, "Contains sphinganine backbone with a fatty acyl group attached to N"