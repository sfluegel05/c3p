"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA is a fatty acyl-CoA that results from the formal condensation
    of the thiol group of coenzyme A with the carboxy group of any long-chain (C13 to C22) fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the thioester functional group pattern (S-C(=O))
    thioester_pattern = Chem.MolFromSmarts("[S][C](=O)")
    matches = mol.GetSubstructMatches(thioester_pattern)
    
    if not matches:
        return False, "No thioester functional group found"
    
    # Assume only one thioester group is present
    s_idx, c_idx = matches[0]
    s_atom = mol.GetAtomWithIdx(s_idx)
    c_atom = mol.GetAtomWithIdx(c_idx)
    
    # Get the bond between sulfur and carbonyl carbon
    bond = mol.GetBondBetweenAtoms(s_idx, c_idx)
    if bond is None:
        return False, "Thioester bond not found"
    
    bond_idx = bond.GetIdx()
    
    # Break the bond to obtain fragments
    fragmented_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    
    # Get the fragments as individual molecules
    fragments = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
    
    # Identify the acyl chain fragment (the one without sulfur)
    acyl_chain = None
    for frag in fragments:
        has_sulfur = any(atom.GetAtomicNum() == 16 for atom in frag.GetAtoms())
        if not has_sulfur:
            acyl_chain = frag
            break
    
    if acyl_chain is None:
        return False, "Acyl chain fragment could not be identified"
    
    # Count the number of carbon atoms in the acyl chain
    carbon_count = sum(1 for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check if the carbon count is between 13 and 22 inclusive
    if 13 <= carbon_count <= 22:
        return True, f"Acyl chain has {carbon_count} carbons, within the long-chain range"
    else:
        return False, f"Acyl chain has {carbon_count} carbons, not within the long-chain range (13-22)"