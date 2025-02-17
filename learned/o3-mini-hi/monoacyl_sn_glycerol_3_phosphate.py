"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate
Definition:
  An sn-glycero-3-phosphate compound having a single unspecified acyl group at either position 1 or position 2.
  
This new approach adds an explicit check on the phosphate group:
  • There must be exactly one phosphorus.
  • The phosphate group must appear as P(=O)(O)(O), and exactly one of its oxygen neighbors must be bound to a carbon.
This helps to reject diesterified phosphates present in e.g. phosphatidylcholines.
Additionally, we count acyl ester groups (excluding those attached to P) and verify the presence of a glycerol-phosphate backbone.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if the molecule given by the SMILES string is a monoacyl-sn-glycerol 3-phosphate.
    
    A monoacyl-sn-glycerol 3-phosphate is defined as an sn-glycero-3-phosphate compound having 
    exactly one acyl group (esterified fatty acid) attached at either position 1 or 2. In addition,
    the phosphate must be monoesterified—that is, it should connect only via one O-C bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Check for exactly one phosphorus atom.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Molecule must contain exactly one phosphorus atom (found {len(phosphorus_atoms)})"
    
    # 3. Check that a phosphate group is present using a SMARTS for P(=O)(O)(O)
    phosphate_smarts = "P(=O)(O)(O)"
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # 4. For the unique phosphorus, ensure it is monoesterified.
    #    Look at all oxygen neighbors of the phosphorus
    p_atom = phosphorus_atoms[0]
    oxy_neighbors = [nbr for nbr in p_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    # Count how many of these oxygens are attached to a carbon (i.e. part of a P–O–C linkage)
    carbon_link_count = 0
    for o in oxy_neighbors:
        # Check neighbors of o (excluding the phosphorus) for a carbon.
        for nbr in o.GetNeighbors():
            if nbr.GetIdx() != p_atom.GetIdx() and nbr.GetAtomicNum() == 6:
                carbon_link_count += 1
                break
    if carbon_link_count != 1:
        return False, f"Phosphate group must be monoesterified (found {carbon_link_count} O-C link(s) on phosphorus)"
    
    # 5. Identify the acyl ester group.
    #    We look for the typical acyl ester fragment: a carbonyl (C(=O)) connected via oxygen to a carbon.
    #    Use SMARTS "[CX3](=O)O[#6]" and ignore cases where the oxygen is attached to a phosphorus.
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)O[#6]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    acyl_count = 0
    for match in acyl_matches:
        # match is a tuple of atom indices: (carbonyl carbon, ester oxygen, alkyl carbon)
        o_atom = mol.GetAtomWithIdx(match[1])
        # Exclude these groups if the oxygen is attached to any phosphorus.
        if any(neighbor.GetAtomicNum() == 15 for neighbor in o_atom.GetNeighbors()):
            continue
        acyl_count += 1
    if acyl_count != 1:
        return False, f"Found {acyl_count} acyl ester group(s); expected exactly one for a monoacyl-sn-glycerol 3-phosphate"
    
    # 6. Check for a glycerol 3-phosphate backbone.
    #    Here we require a simplified pattern for the glycerol moiety attached to the phosphate.
    #    The pattern "OC(CO)COP(=O)(O)O" (ignoring stereochemistry) is chosen to match the expected backbone.
    glycerol_phosphate_smarts = "OC(CO)COP(=O)(O)O"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol 3-phosphate backbone not found"
    
    # If all tests pass, the molecule meets the criteria.
    return True, "Molecule has a glycerol 3-phosphate backbone with exactly one acyl ester group, and the phosphate is monoesterified"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Try one of the provided SMILES strings.
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O"  # 1-icosanoyl-sn-glycero-3-phosphate
    result, reason = is_monoacyl_sn_glycerol_3_phosphate(test_smiles)
    print(result, reason)