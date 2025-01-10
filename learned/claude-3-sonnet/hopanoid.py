"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    Hopanoids are pentacyclic triterpenoids based on the hopane skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic hopane skeleton SMARTS - more flexible pattern
    # Accounts for both saturated and unsaturated versions
    hopane_core = Chem.MolFromSmarts("[C]1[C,c]2[C,c][C,c][C,c]3[C,c]([C,c]2[C,c][C,c][C,c]1)[C,c][C,c][C,c]4[C,c]3[C,c][C,c][C,c]4")
    
    if not mol.HasSubstructMatch(hopane_core):
        return False, "No hopane pentacyclic core found"

    # Count rings - hopanoids should have at least 5 rings
    ri = mol.GetRingInfo()
    if ri.NumRings() < 5:
        return False, "Insufficient number of rings"

    # Count carbons - hopanoids typically have 27-35 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27 or c_count > 35:
        return False, f"Carbon count ({c_count}) outside typical hopanoid range (27-35)"

    # Check for characteristic methyl groups
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_count < 6:
        return False, f"Insufficient methyl groups ({methyl_count}), hopanoids typically have 6+ methyl groups"

    # Verify basic molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350 or mol_wt > 1000:
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical hopanoid range"

    return True, "Contains hopane pentacyclic core with characteristic features"