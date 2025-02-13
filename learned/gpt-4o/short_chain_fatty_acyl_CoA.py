"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

    A short-chain fatty acyl-CoA is a fatty acyl-CoA formed by the condensation
    of the thiol group of coenzyme A with a short-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive pattern for Coenzyme A
    coenzymeA_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)CCSCCC(=O)NCCSC(=S)CC(O)CN1C=NC2=C1N=CN=C2N")
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Define the thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Ensure short-chain fatty acid component (2-6 carbons)
    fatty_acid_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and bond.GetBeginAtom().GetSymbol() == 'S' and bond.GetEndAtom().GetSymbol() == 'C':
            # Get the carbon chain from the sulfur atom
            carbon_atom = bond.GetEndAtom() if bond.GetBeginAtom().GetSymbol() == 'S' else bond.GetBeginAtom()
            visited = set()
            carbon_count = 0
            to_visit = [carbon_atom]
            while to_visit:
                atom = to_visit.pop()
                if atom.GetIdx() not in visited and atom.GetSymbol() == 'C':
                    visited.add(atom.GetIdx())
                    carbon_count += 1
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in visited:
                            to_visit.append(neighbor)
            if 2 <= carbon_count <= 6:
                fatty_acid_found = True
                break

    if not fatty_acid_found:
        return False, "No proper short-chain fatty acyl component found"
    
    # If all checks pass, classify as short-chain fatty acyl-CoA
    return True, "Contains a CoA moiety joined by thioester bond to a short-chain fatty acid"