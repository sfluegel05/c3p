"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: CHEBI:73048 O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is an O-acylcarnitine with the carnitine component in the L-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carnitine backbone with L-configuration
    # More flexible pattern to match the core structure
    carnitine_pattern = Chem.MolFromSmarts("[N+X4](C)(C)C[C@H](CC([O-])=O)OC(=O)")
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "No L-carnitine backbone found"

    # Check for the presence of a quaternary ammonium group (N+)
    ammonium_pattern = Chem.MolFromSmarts("[N+X4]")
    if not mol.HasSubstructMatch(ammonium_pattern):
        return False, "No quaternary ammonium group found"

    # Check for the carboxylate group (COO-)
    carboxylate_pattern = Chem.MolFromSmarts("CC([O-])=O")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for the ester bond (O-acyl group)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester bond found"

    # Check for the L-configuration by ensuring the chiral center is correctly configured
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral center found"
    
    # Ensure the chiral center is in the L-configuration
    # L-configuration corresponds to 'R' in RDKit for carnitine derivatives
    for center in chiral_centers:
        if center[1] != 'R':  # Changed from 'S' to 'R'
            return False, "Chiral center not in L-configuration"

    # Additional check to reduce false positives: ensure the molecule has a reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 1500:  # Expanded weight range
        return False, "Molecular weight out of range for O-acyl-L-carnitine"

    # Additional check: count the number of carbons to ensure it's a reasonable acyl chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7 or c_count > 50:  # Reasonable range for acyl chains
        return False, "Number of carbons out of range for O-acyl-L-carnitine"

    return True, "Contains L-carnitine backbone with O-acyl group, quaternary ammonium, and carboxylate group"