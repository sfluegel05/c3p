"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:35748 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (-CONH- or -CONR-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H0,H1]")
    amide_match = mol.GetSubstructMatch(amide_pattern)
    if not amide_match:
        return False, "No amide group found"

    # Get indices of atoms in the amide group
    amide_atom_indices = set(amide_match)

    # Count carbons not in the amide group
    c_count = sum(1 for atom in mol.GetAtoms() 
                 if atom.GetAtomicNum() == 6 and 
                 atom.GetIdx() not in amide_atom_indices)
    
    if c_count < 6:
        return False, f"Not enough carbons ({c_count}) for a fatty chain"

    # Check molecular weight - fatty amides typically >150 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a fatty amide"

    # Check for a continuous hydrocarbon chain of at least 6 carbons
    # This pattern allows for branches and unsaturations
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No suitable hydrocarbon chain found"

    # Additional check: at least 50% of carbons should be in chains
    # (excluding amide group carbons)
    total_c = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count/total_c < 0.5:
        return False, "Too many carbons in non-chain positions"

    # Check that the amide group is connected to the hydrocarbon chain
    # Get the carbon in the amide group
    amide_carbon = [atom for atom in amide_match if atom.GetAtomicNum() == 6 and atom.GetDegree() == 3][0]
    
    # Check if the amide carbon is connected to at least one carbon chain
    has_chain = any(neighbor.GetAtomicNum() == 6 and 
                   neighbor.GetIdx() not in amide_atom_indices 
                   for neighbor in amide_carbon.GetNeighbors())
    
    if not has_chain:
        return False, "Amide group not connected to a hydrocarbon chain"

    return True, "Contains an amide group connected to a long hydrocarbon chain"