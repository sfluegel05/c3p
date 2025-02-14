"""
Classifies: CHEBI:29348 fatty amide
"""
"""
Classifies: CHEBI:33670 fatty amide
A fatty amide is defined as a monocarboxylic acid amide derived from a fatty acid.
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
    amide_pattern = Chem.MolFromSmarts("[NX3](=[OX1])([#6])[#6]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"
    
    # Look for fatty acid chain (long carbon chain)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = []
    for amide_idx in amide_matches:
        amide_atom = mol.GetAtomWithIdx(amide_idx)
        for bond in amide_atom.GetBonds():
            neighbor_idx = bond.GetOtherAtomIdx(amide_idx)
            if neighbor_idx != -1 and mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 6:
                fatty_acid_match = mol.GetAtomWithIdx(neighbor_idx).IsInRingSize(0) and \
                                   mol.HasSubstructMatch(fatty_acid_pattern, AtomIds=[neighbor_idx])
                if fatty_acid_match:
                    fatty_acid_matches.append((amide_idx, neighbor_idx))
                    break
    
    if not fatty_acid_matches:
        return False, "No fatty acid chain attached to amide group"
    
    # Check molecular weight and rotatable bonds for fatty acid chain
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if mol_wt < 200 or n_rotatable < 5:
        return False, "Fatty acid chain is too short"
    
    # Additional checks to exclude false positives
    unwanted_patterns = ['c1ccccc1', 'O=C(O)C', 'C(=O)O', 'N#C', 'S(=O)(=O)', 'P']
    for pattern in unwanted_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains unwanted substructure"
    
    # Check the length of the fatty acid chain
    for amide_idx, chain_start_idx in fatty_acid_matches:
        chain_length = 0
        curr_atom_idx = chain_start_idx
        while mol.GetAtomWithIdx(curr_atom_idx).GetAtomicNum() == 6:
            chain_length += 1
            for bond in mol.GetAtomWithIdx(curr_atom_idx).GetBonds():
                neighbor_idx = bond.GetOtherAtomIdx(curr_atom_idx)
                if neighbor_idx != -1 and mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 6:
                    curr_atom_idx = neighbor_idx
                    break
        if chain_length < 8:
            return False, "Fatty acid chain is too short"
    
    return True, "Contains amide group with a fatty acid chain attached"