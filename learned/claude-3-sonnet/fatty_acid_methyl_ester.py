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
    
    # Count number of methyl ester groups
    num_methyl_esters = len(methyl_ester_matches)
    if num_methyl_esters > 2:
        return False, f"Too many methyl ester groups ({num_methyl_esters})"
    
    # Get atoms in methyl ester group
    ester_atoms = set()
    for match in methyl_ester_matches:
        ester_atoms.update(match)
    
    # Count rings and check their nature
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    # Limit total number of rings
    if num_rings > 3:
        return False, "Too many ring systems for a fatty acid"
    
    # Check for aromatic rings
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings > 0:
        return False, "Fatty acid methyl esters should not contain aromatic rings"
    
    # Count carbons and check chain
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    non_ester_carbons = total_carbons - (2 * num_methyl_esters)  # Subtract carbons from methyl ester groups
    
    if non_ester_carbons < 3:
        return False, "Carbon chain too short for fatty acid"
    
    if non_ester_carbons > 30:
        return False, "Carbon chain too long for typical fatty acid"
    
    # Check for predominantly aliphatic character
    num_sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if num_sp3_carbons < (total_carbons * 0.5):
        return False, "Not predominantly aliphatic"
    
    # Special case for dimethyl esters (like dimethyl sebacate)
    if num_methyl_esters == 2:
        # Check for reasonable chain length between esters
        if non_ester_carbons < 6:
            return False, "Chain too short for dicarboxylic acid"
        return True, "Valid dimethyl ester of dicarboxylic acid"
    
    # Check for continuous carbon chain
    longest_chain = rdMolDescriptors.CalcNumBonds(mol)
    if longest_chain < 5:  # Minimum chain length including the ester group
        return False, "No proper fatty acid chain"
    
    # Count heteroatoms (excluding the ester oxygens)
    num_N = rdMolDescriptors.CalcNumLipinskiHBA(mol) - (2 * num_methyl_esters)
    num_S = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[#16]")))
    if num_N > 1 or num_S > 1:
        return False, "Too many heteroatoms for a fatty acid"
    
    return True, "Contains methyl ester group with appropriate fatty acid chain"