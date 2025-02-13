"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:35506 fatty acyl-CoA
An acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any fatty acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for coenzyme A backbone pattern
    coa_pattern = Chem.MolFromSmarts("[C@H]1([C@@H](O[P@@](O)(=O)O[P@](O)(=O)OC[C@H]2[C@@H](N3C=NC4=C3N=CN=C4N)[C@H](O)[C@@H](O)[C@H]2O[P](O)(O)=O)O)[C@H](O)[C@H](O)[C@@H]1OP(O)(O)=O")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No coenzyme A backbone found"

    # Look for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Look for fatty acid chain (long carbon chain attached to thioester)
    fatty_acid_patterns = [
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2]"),  # Saturated chains
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]"),
        Chem.MolFromSmarts("[CH3]C(C)(C)"),  # Branched chains
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(C)(C)"),
        Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2]C(C)(C)[CH2]"),
        Chem.MolFromSmarts("[CH2][CH2]=C"),  # Unsaturated chains
        Chem.MolFromSmarts("[CH2][CH2]=C[CH2][CH2]=C"),
        Chem.MolFromSmarts("[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C"),
        Chem.MolFromSmarts("[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C"),
        Chem.MolFromSmarts("[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C"),
        Chem.MolFromSmarts("[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C[CH2][CH2]=C"),
    ]
    fatty_acid_matches = []
    for pattern in fatty_acid_patterns:
        fatty_acid_matches.extend(mol.GetSubstructMatches(pattern))
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"

    # Check for additional functional groups (e.g., hydroxy, keto, etc.)
    additional_groups_pattern = Chem.MolFromSmarts("[OH,O=C]")
    additional_groups_matches = mol.GetSubstructMatches(additional_groups_pattern)

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - fatty acyl-CoAs typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for fatty acyl-CoA"

    # Check for additional structural rules
    if additional_groups_matches:
        return True, "Contains coenzyme A backbone, thioester linkage, fatty acid chain, and additional functional groups"
    else:
        return True, "Contains coenzyme A backbone, thioester linkage, and fatty acid chain"