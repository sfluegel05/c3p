"""
Classifies: CHEBI:27093 tricarboxylic acid
"""
"""
Classifies: tricarboxylic acid
Definition: An oxoacid containing three carboxy groups
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tricarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid based on its SMILES string.
    A tricarboxylic acid must contain exactly three carboxylic acid (-COOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except:
        return False, "Molecule failed sanitization"

    # Pattern that matches both protonated and deprotonated carboxylic acids
    carboxyl_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2]')
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    
    if len(matches) < 3:
        return False, f"Found only {len(matches)} carboxyl groups, need exactly 3"
    
    # Get all matching carbons and oxygens for further verification
    carboxyl_carbons = set(match[0] for match in matches)
    carboxyl_oxygens = set(match[2] for match in matches)
    
    # Count actual carboxylic acid groups by verifying hydrogens or charge
    true_acid_count = 0
    for carbon_idx in carboxyl_carbons:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        oxygen_neighbors = [n for n in carbon.GetNeighbors() 
                          if n.GetAtomicNum() == 8 and n.GetIdx() in carboxyl_oxygens]
        
        for oxygen in oxygen_neighbors:
            # Check if oxygen is part of carboxylic acid
            if (oxygen.GetTotalNumHs() > 0 or oxygen.GetFormalCharge() == -1) and \
               not any(n.GetAtomicNum() == 6 and n.GetIdx() != carbon_idx for n in oxygen.GetNeighbors()):
                true_acid_count += 1
                break
    
    if true_acid_count != 3:
        return False, f"Found {true_acid_count} true carboxylic acid groups, need exactly 3"

    # Check for potentially interfering groups
    ester_pattern = Chem.MolFromSmarts('[CX3](=[OX1])[OX2][CX4]')
    ester_matches = set(match[0] for match in mol.GetSubstructMatches(ester_pattern))
    
    # Only count esters that share carbons with carboxylic acids
    interfering_esters = ester_matches.intersection(carboxyl_carbons)
    if interfering_esters:
        return False, "Contains ester groups that share carbons with carboxylic acids"
    
    # Check molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:  # Minimum reasonable weight for tricarboxylic acid
        return False, "Molecular weight too low for tricarboxylic acid"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for tricarboxylic acid"
    if o_count < 6:
        return False, "Too few oxygens for tricarboxylic acid"

    # Verify carboxylic groups are independent
    for i, c1 in enumerate(carboxyl_carbons):
        for c2 in list(carboxyl_carbons)[i+1:]:
            path = Chem.GetShortestPath(mol, c1, c2)
            if len(path) < 3:  # Carbons directly connected or sharing atom
                return False, "Carboxylic acid groups are not independent"

    return True, "Contains exactly three independent carboxylic acid groups"