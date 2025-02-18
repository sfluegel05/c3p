"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: CHEBI:37237 N-acylphytosphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phytosphingosine backbone pattern
    # (1) Long carbon chain (>= 16C)
    # (2) OH groups at C1, C3, C4 
    # (3) NH group at C2
    phytosphingosine_pattern = Chem.MolFromSmarts("[C;H3][C;H2][C;H2][C;H2]([NH][C;H2][C;H1]([C;H3])[C;H2][C;H2][C;H2]([C;H3])[C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H3])[OH][C;H2][OH][C;H2][OH]")
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
        return False, "No phytosphingosine backbone found"
        
    # Look for acyl group attached to nitrogen (-N-C(=O)-)
    acyl_pattern = Chem.MolFromSmarts("[N;H1][C;X3](=[O;X1])")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) == 0:
        return False, "No acyl group attached to nitrogen"

    # Check for fatty acid chain (long carbon chain attached to acyl)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern, acyl_matches[0][0])
    if len(fatty_acid_matches) == 0:
        return False, "No fatty acid chain attached to acyl group"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for N-acylphytosphingosine"

    return True, "Contains phytosphingosine backbone with fatty acyl group attached to nitrogen"