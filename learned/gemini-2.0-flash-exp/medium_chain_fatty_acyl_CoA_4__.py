"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the core CoA structure (more generalized)
    # This pattern focuses on key functional groups: pantetheine, pyrophosphate, ribose-adenine
    coa_pattern = Chem.MolFromSmarts("NCC(=O)CCNC(=O)[C](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C]1O[C]([C](O)[C](O)1)n1cnc2c(N)ncnc12") # Removed specific stereochemistry
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA substructure not found"
    
    # Check for thioester bond (-C(=O)-S-) connected to CoA substructure
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Find the acyl chain attached to the thioester using a more specific SMARTS
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)S-[CX4]")
    acyl_chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if not acyl_chain_matches:
        return False, "No acyl chain attached to thioester"
    
    # Count the carbons in acyl chain (excluding the carbonyl carbon)
    acyl_chain_start_atom_idx = acyl_chain_matches[0][2]  # Third atom in the match (carbon attached to sulfur)
    
    # Function to count carbons in the acyl chain
    def count_chain_carbons(mol, start_atom_idx):
      
        carbons = 0
        
        queue = [start_atom_idx]
        visited = set()
        
        while queue:
          current_atom_idx = queue.pop(0)
          
          if current_atom_idx in visited:
            continue
          visited.add(current_atom_idx)

          current_atom = mol.GetAtomWithIdx(current_atom_idx)
          
          if current_atom.GetAtomicNum() == 6:
            carbons += 1
            for neighbor in current_atom.GetNeighbors():
               if neighbor.GetAtomicNum() == 6:
                queue.append(neighbor.GetIdx())
        return carbons
            
    chain_carbons = count_chain_carbons(mol,acyl_chain_start_atom_idx)
    
    
    if  chain_carbons < 6 or chain_carbons > 12:
           return False, f"Acyl chain has {chain_carbons} carbons, must be between 6 and 12"
           
    # Verify the total charge of -4 using the formal charges of the atoms
    total_charge = 0
    for atom in mol.GetAtoms():
      total_charge += atom.GetFormalCharge()

    if total_charge != -4:
        return False, f"Molecule has {total_charge} charge instead of -4"

    return True, "Matches medium-chain fatty acyl-CoA(4-) criteria"