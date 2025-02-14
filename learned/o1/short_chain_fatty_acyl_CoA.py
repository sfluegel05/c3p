"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60940 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA results from the condensation of the thiol group of coenzyme A
    with the carboxy group of any short-chain fatty acid (2 to 5 carbons long).
    
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
    
    # Define a general SMARTS pattern for Coenzyme A moiety focusing on adenine ring and phosphates
    # This pattern matches the adenine ring attached to ribose with phosphates
    coA_smarts_str = "n1cnc2c(ncnc12)[C@H]3O[C@H](COP(O)(=O)O)C(O)[C@@H]3OP(O)(O)=O"
    coA_smarts = Chem.MolFromSmarts(coA_smarts_str)
    if coA_smarts is None:
        return False, "Invalid CoA SMARTS pattern"
    
    # Check for Coenzyme A substructure without considering stereochemistry
    if not mol.HasSubstructMatch(coA_smarts, useChirality=False):
        return False, "Coenzyme A moiety not found"
    
    # Define SMARTS pattern for thioester linkage (C(=O)-S)
    thioester_smarts_str = "C(=O)S"
    thioester_smarts = Chem.MolFromSmarts(thioester_smarts_str)
    if thioester_smarts is None:
        return False, "Invalid thioester SMARTS pattern"
    
    # Find the thioester linkage
    thioester_matches = mol.GetSubstructMatches(thioester_smarts)
    if not thioester_matches:
        return False, "Thioester linkage not found"
    
    # For each thioester linkage, attempt to identify the acyl chain and CoA moiety
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbon of C=O
        sulfur_idx = match[2]      # Sulfur atom index
        
        # Get atoms connected to the carbonyl carbon (should be O and alpha carbon of acyl chain)
        carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)
        neighbors_c = [a.GetIdx() for a in carbonyl_c.GetNeighbors() if a.GetIdx() != sulfur_idx]
        
        # Get atoms connected to the sulfur atom (should be connected to CoA moiety)
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        neighbors_s = [a.GetIdx() for a in sulfur_atom.GetNeighbors() if a.GetIdx() != carbonyl_c_idx]
        
        # Create fragments by breaking the C-S bond of the thioester linkage
        bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, sulfur_idx)
        mol_frag = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=False)
        frags = Chem.GetMolFrags(mol_frag, asMols=True, sanitizeFrags=True)
        
        # Identify acyl chain and CoA fragments
        acyl_chain = None
        coa_fragment = None
        for frag in frags:
            # If fragment contains the carbonyl carbon, it's the acyl chain
            if any(atom.GetAtomicNum() == 6 and atom.GetIdx() == carbonyl_c_idx for atom in frag.GetAtoms()):
                acyl_chain = frag
            # If fragment contains the sulfur connected to CoA, it's the CoA fragment
            elif any(atom.GetAtomicNum() == 16 and atom.GetIdx() == sulfur_idx for atom in frag.GetAtoms()):
                coa_fragment = frag
        
        if acyl_chain is None or coa_fragment is None:
            continue
        
        # Check if the CoA fragment contains the CoA SMARTS pattern
        if not coa_fragment.HasSubstructMatch(coA_smarts, useChirality=False):
            continue  # Skip if CoA pattern not found in this fragment
        
        # Count the number of carbons in the acyl chain excluding the carbonyl carbon
        acyl_carbons = [atom for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6]
        num_carbons = len(acyl_carbons) - 1  # Exclude the carbonyl carbon
        
        if 1 <= num_carbons <= 4:
            return True, f"Found short-chain fatty acyl-CoA with acyl chain length {num_carbons + 1} carbons"
        else:
            return False, f"Acyl chain length is {num_carbons + 1}, not short-chain"

    return False, "Acyl chain not identified or not short-chain"

# Add metadata
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:60940',
        'name': 'short-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any short-chain fatty acid.',
        'parents': ['CHEBI:37554', 'CHEBI:57288']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}