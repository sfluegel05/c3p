"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
            
        # Get basic molecular properties
        c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
        rings = rdMolDescriptors.CalcNumRings(mol)
        
        # Check carbon count - should be around 20 (allowing some variation due to modifications)
        if c_count < 15 or c_count > 40:
            return False, f"Carbon count ({c_count}) outside typical range for diterpenoids (15-40)"
            
        # Diterpenoids typically have multiple rings
        if rings < 2:
            return False, f"Too few rings ({rings}) for a diterpenoid"
            
        # Check for typical terpenoid characteristics
        # Look for common structural features in diterpenoids
        
        # Check for branched carbon chains
        branched_carbon = Chem.MolFromSmarts("[CH1,CH0](-[CH3])(-[CH2,CH3])")
        if not mol.HasSubstructMatch(branched_carbon):
            return False, "Missing characteristic branched carbon structure"
            
        # Look for ring systems common in diterpenoids
        decalin_like = Chem.MolFromSmarts("C1CCC2CCCCC2C1")  # Basic decalin-like structure
        
        # Check molecular weight - most diterpenoids are between 250-500 Da
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        if mol_wt < 250 or mol_wt > 900:
            return False, f"Molecular weight ({mol_wt:.1f}) outside typical range for diterpenoids"
            
        # Count sp3 and sp2 carbons to check for degree of unsaturation
        sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
        sp2_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3]")))
        
        # Most diterpenoids have a mix of sp3 and sp2 carbons
        if sp3_carbons < 5:
            return False, "Too few sp3 carbons for a diterpenoid"
            
        # Check for presence of common functional groups in diterpenoids
        hydroxy = mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))
        ketone = mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[C,c]"))
        ester = mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O][C,c]"))
        double_bond = mol.HasSubstructMatch(Chem.MolFromSmarts("[C]=[C]"))
        
        # Count number of these features
        functional_groups = sum([hydroxy, ketone, ester, double_bond])
        
        # Most diterpenoids have at least one of these functional groups
        if functional_groups == 0:
            return False, "Missing typical diterpenoid functional groups"
            
        # Calculate degree of unsaturation (rings + double bonds)
        double_bonds = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
        saturation_degree = rings + double_bonds
        
        if saturation_degree < 2:
            return False, "Insufficient rings and double bonds for a diterpenoid"
            
        # If we've made it here, check for overall structural characteristics
        if (c_count >= 15 and c_count <= 40 and 
            rings >= 2 and 
            functional_groups > 0 and
            250 <= mol_wt <= 900):
            return True, "Matches diterpenoid characteristics: carbon count, ring systems, functional groups, and molecular weight"
            
        return False, "Does not match overall diterpenoid characteristics"
        
    except Exception as e:
        return None, f"Error in analysis: {str(e)}"