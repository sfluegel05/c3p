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
    
    # Get the carbonyl carbon atom index
    carbonyl_carbon_idx = methyl_ester_matches[0][0]
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_idx)
    
    # Check that there's at least one carbon attached to the carbonyl carbon
    # that's not part of the methyl ester group
    non_ester_neighbors = [neighbor for neighbor in carbonyl_carbon.GetNeighbors() 
                         if neighbor.GetAtomicNum() == 6]
    if not non_ester_neighbors:
        return False, "No carbon chain attached to ester group"
    
    # Count carbons not in the methyl ester group
    methyl_ester_atoms = set([atom_idx for match in methyl_ester_matches for atom_idx in match])
    non_ester_carbons = sum(1 for atom in mol.GetAtoms() 
                           if atom.GetAtomicNum() == 6 and atom.GetIdx() not in methyl_ester_atoms)
    
    if non_ester_carbons < 3:
        return False, "Carbon chain too short for fatty acid"
    
    # Check for aromatic rings - fatty acid methyl esters are typically aliphatic
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings > 0:
        return False, "Fatty acid methyl esters should not contain aromatic rings"
    
    # Allow some rings (can have epoxy groups etc)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 5:  # Increased threshold to allow more complex structures
        return False, "Too many ring systems"
        
    # Check that the carbon attached to ester is sp2 (carbonyl)
    if carbonyl_carbon.GetHybridization() != Chem.HybridizationType.SP2:
        return False, "Ester carbon must be sp2 hybridized"
    
    # Count heteroatoms (excluding O from ester group)
    num_N = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if num_N > 0:
        return False, "Should not contain nitrogen atoms"
        
    # Allow various modifications (hydroxy, epoxy, peroxy, halogens, etc)
    # No restrictions on halogens or other modifications as they can be present
    
    return True, "Contains single methyl ester group with appropriate fatty acid chain"