"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids include sphinganine, its homologs and stereoisomers, 
    and hydroxyl or unsaturated derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the 2-amino-1,3-diol backbone characteristic of sphingoids
    amino_diol_pattern = Chem.MolFromSmarts("N[C@@H](CO)C(O)")
    if not mol.HasSubstructMatch(amino_diol_pattern):
        return False, "Missing core 2-amino-1,3-diol moiety"
    
    # Look for long hydrocarbon chain (at least 14 carbons)
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No sufficiently long hydrocarbon chain found"
    
    # Check for either unsaturation or additional hydroxyl groups
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    
    has_double_bond = mol.HasSubstructMatch(double_bond_pattern)
    has_additional_hydroxyl = len(mol.GetSubstructMatches(hydroxyl_pattern)) > 1  # More than one hydroxyl
    if not (has_double_bond or has_additional_hydroxyl):
        return False, "Lacks unsaturation or additional hydroxyl groups, at least one is needed"

    # Check for stereochemistry typical of sphingoids
    stereo_moieties = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(stereo_moieties) == 0:
        return False, "No stereochemistry typical of sphingoids found"

    return True, "SMILES corresponds to a recognized sphingoid structure"