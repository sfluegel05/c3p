"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoids
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are characterized by specific structural features such as:
    - The dibenzopyran ring system or variants in phytocannabinoids (e.g., THC, CBD).
    - Long-chain fatty acid derivatives linked to ethanolamine or glycerol (endocannabinoids).
    - Synthetic cannabinoids with specific core structures and side chains.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Normalize molecule (e.g., handle tautomers)
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    
    # Define patterns for phytocannabinoids (e.g., THC, CBD)
    # These patterns are designed to match the core structures of THC and CBD with variability
    thc_cbd_pattern = Chem.MolFromSmarts("""
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]-2=[#6]-[#6]-[#6]-[#6]-2
    """)
    
    # Alternative pattern for cannabinoids with alkyl side chain
    cannabinoid_core_pattern = Chem.MolFromSmarts("""
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-[#6]2[#6][#6][#6][#6][#6]2
    """)
    
    # Check for phytocannabinoid features
    if mol.HasSubstructMatch(thc_cbd_pattern) or mol.HasSubstructMatch(cannabinoid_core_pattern):
        return True, "Contains core structure characteristic of phytocannabinoids like THC and CBD"
    
    # Define patterns for endocannabinoids (N-acylethanolamines and monoacylglycerols)
    n_acylethanolamine_pattern = Chem.MolFromSmarts('C(=O)NCC[OX2H]')  # Amide linked to ethanolamine
    monoacylglycerol_pattern = Chem.MolFromSmarts('C(=O)OCC(O)CO')     # Ester linked to glycerol
    
    # Check for endocannabinoid features (N-acylethanolamines)
    if mol.HasSubstructMatch(n_acylethanolamine_pattern):
        # Analyze fatty acid chain length without strict unsaturation criteria
        carbonyl_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetTotalValence() == 3 and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetBeginAtomIdx() == atom.GetIdx() for bond in atom.GetBonds())]
        for carbon_idx in carbonyl_carbons:
            chain_length = _count_chain_length(mol, carbon_idx)
            if chain_length >= 12:
                return True, "Contains N-acylethanolamine structure characteristic of endocannabinoids"
        # If no sufficient chain found
        return False, "Amide linked to ethanolamine found but fatty acid chain too short"
    
    # Check for endocannabinoid features (monoacylglycerols)
    if mol.HasSubstructMatch(monoacylglycerol_pattern):
        # Analyze fatty acid chain length without strict unsaturation criteria
        carbonyl_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetTotalValence() == 3 and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetBeginAtomIdx() == atom.GetIdx() for bond in atom.GetBonds())]
        for carbon_idx in carbonyl_carbons:
            chain_length = _count_chain_length(mol, carbon_idx)
            if chain_length >= 12:
                return True, "Contains monoacylglycerol structure characteristic of endocannabinoids"
        # If no sufficient chain found
        return False, "Ester linked to glycerol found but fatty acid chain too short"
    
    # Define refined pattern for synthetic cannabinoids (e.g., specific indole derivatives)
    synthetic_cannabinoid_patterns = [
        Chem.MolFromSmarts('c1nccc2c1cccc2'),  # Indole core
        Chem.MolFromSmarts('c1cncc2c1cccc2'),  # Indazole core
        Chem.MolFromSmarts('c1cc(cnc1)c2ccc(cc2)C(=O)'),  # Indole with acyl side chain
        Chem.MolFromSmarts('c1ccc(cc1)C(=O)N2CCCCC2'),  # Naphthoylindole
        Chem.MolFromSmarts('c1ccc2c(c1)ccc(c2)C(=O)N3CCCCC3')  # More specific naphthoylindole
    ]
    
    # Check for synthetic cannabinoid features
    for pattern in synthetic_cannabinoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains core structure characteristic of synthetic cannabinoids"
    
    # Additional check for cannabigerol-like structures
    cannabigerol_pattern = Chem.MolFromSmarts('CCCCCc1cc(O)c(cc1O)CCC=C(C)C')
    if mol.HasSubstructMatch(cannabigerol_pattern):
        return True, "Contains structure characteristic of cannabigerol-type cannabinoids"
    
    return False, "No characteristic cannabinoid structural features found"

def _count_chain_length(mol, start_idx):
    """
    Counts the number of carbons in the fatty acid chain starting from the given atom index.
    Args:
        mol: RDKit Mol object
        start_idx: Index of the starting carbon atom
    Returns:
        int: Number of carbons in the chain
    """
    visited = set()
    to_visit = [start_idx]
    chain_length = 0
    while to_visit:
        current_idx = to_visit.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)
        atom = mol.GetAtomWithIdx(current_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain_length += 1
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                    to_visit.append(neighbor_idx)
    return chain_length