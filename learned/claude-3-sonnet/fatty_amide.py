"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: fatty amide compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Look for amide group
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    if not amide_matches:
        return False, "No amide group found"
    
    # Get number of amide groups
    num_amides = len(amide_matches)
    
    # Check carbon chain length
    carbon_chain = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(carbon_chain)
    
    if not chain_matches:
        return False, "No significant carbon chain found"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Carbon chain too short for fatty amide"
        
    # Count rotatable bonds to check chain flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Additional checks for specific cases that should be excluded
    
    # Check for cyclic peptides (multiple amides in cycle)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        amide_in_ring = 0
        for amide_match in amide_matches:
            if amide_match[0] in ring:
                amide_in_ring += 1
        if amide_in_ring > 1:
            return False, "Appears to be a cyclic peptide"
    
    # Check for peptides (multiple amides in sequence)
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][NX3][CX3](=[OX1])")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Appears to be a peptide"
        
    # Positive classification criteria
    if num_amides == 1:
        if c_count >= 4 and n_rotatable >= 2:
            return True, "Contains single amide group with appropriate carbon chain"
    else:
        # For molecules with multiple amides, check if they're separated by long chains
        if c_count/num_amides >= 4 and n_rotatable >= num_amides*2:
            return True, "Contains well-separated amide groups with appropriate carbon chains"
            
    return False, "Does not match fatty amide criteria"