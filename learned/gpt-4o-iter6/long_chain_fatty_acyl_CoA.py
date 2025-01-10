"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA includes a long-chain fatty acid esterified with the thiol group of coenzyme A.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the Coenzyme A motif more loosely, capturing the adenine and ribose phosphate components
    ribose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@@H]([C@H](O)[C@@H]([C@@H]1O)O)OP(O)(O)=O")
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    
    if not (mol.HasSubstructMatch(ribose_pattern) and
            mol.HasSubstructMatch(adenine_pattern) and
            mol.HasSubstructMatch(phosphate_pattern)):
        return False, "No Coenzyme A moiety found"

    # Define thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"
    
    # Attempt to find the carbon chain length from the thioester linkage
    for match in thioester_matches:
        carbon_chain_atoms = set()
        # Start from the carbonyl carbon attached to sulfur
        carbon = mol.GetAtomWithIdx(match[0]) 
        visited_atoms = set()
        queue = [carbon]
        
        while queue:
            current_atom = queue.pop()
            if current_atom.GetIdx() not in visited_atoms:
                visited_atoms.add(current_atom.GetIdx())
                # Count contiguous C's moving away from the thioester oxygen
                if current_atom.GetSymbol() == 'C' and current_atom.GetIdx() != match[-1]:
                    carbon_chain_atoms.add(current_atom.GetIdx())
                # Explore unvisited neighbors
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited_atoms:
                        queue.append(neighbor)

        carbon_count = len(carbon_chain_atoms)
        
        if 13 <= carbon_count <= 22:
            return True, f"Valid long-chain fatty acyl-CoA with {carbon_count} C atoms in chain"

    return False, "Carbon chain length not in C13-C22 for long-chain fatty acid"