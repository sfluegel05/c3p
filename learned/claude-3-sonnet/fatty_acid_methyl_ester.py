"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: fatty acid methyl ester
A fatty acid ester that is the carboxylic ester obtained by the formal condensation 
of a fatty acid with methanol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for methyl ester group pattern
    methyl_ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3X4]")
    methyl_ester_matches = mol.GetSubstructMatches(methyl_ester_pattern)
    
    if not methyl_ester_matches:
        return False, "No methyl ester group found"
    
    # Must have exactly one methyl ester group
    if len(methyl_ester_matches) != 1:
        return False, f"Found {len(methyl_ester_matches)} methyl ester groups, need exactly 1"
    
    # Count carbons
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    non_ester_carbons = total_carbons - 2  # Subtract carbons from methyl ester group
    
    if non_ester_carbons < 4:
        return False, "Carbon chain too short for fatty acid"
    
    if non_ester_carbons > 30:
        return False, "Carbon chain too long for typical fatty acid"
        
    # Check for reasonable molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for fatty acid methyl esters"
    
    # Check for aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings > 0:
        return False, "Fatty acid methyl esters should not contain aromatic rings"
    
    # Allow some rings (can have epoxy groups etc)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 4:
        return False, "Too many ring systems"
        
    # Check that the carbon attached to ester is sp2 (carbonyl)
    ester_match = methyl_ester_matches[0]
    carbonyl_carbon = mol.GetAtomWithIdx(ester_match[0])
    if carbonyl_carbon.GetHybridization() != Chem.HybridizationType.SP2:
        return False, "Ester carbon must be sp2 hybridized"
    
    # Check for carbon chain attached to ester
    carbon_chain_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH3X4]~[#6]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No carbon chain attached to ester group"
    
    # Allow various modifications (hydroxy, epoxy, peroxy, etc) but limit halogens
    num_halogens = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[F,Cl,I]")))
    num_Br = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[Br]")))
    if num_halogens > 0 and num_Br == 0:  # Allow Br but not other halogens
        return False, "Contains non-bromine halogens"
    if num_Br > 1:
        return False, "Contains too many bromine atoms"
        
    return True, "Contains single methyl ester group with appropriate fatty acid chain"