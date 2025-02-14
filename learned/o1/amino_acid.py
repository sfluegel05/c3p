"""
Classifies: CHEBI:33709 amino acid
"""
"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is defined as a carboxylic acid containing one or more amino groups.
    It should not be a peptide (no peptide bonds connecting amino acid residues).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    carboxylic_acid_smarts = "[CX3](=O)[O;H1,-]"  # Carboxylic acid (include deprotonated form)
    amino_group_smarts = "[NX3;H2,H1;!$(N-*);!$([N+])]"  # Primary amine (neutral or protonated)
    peptide_bond_smarts = "[NX3][CX3](=O)"  # Peptide bond: N-C(=O)
    
    # Create pattern molecules
    carboxylic_acid_pattern = Chem.MolFromSmarts(carboxylic_acid_smarts)
    amino_group_pattern = Chem.MolFromSmarts(amino_group_smarts)
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    
    # Search for carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    else:
        num_carboxylic_acids = len(carboxylic_acid_matches)
    
    # Search for amino groups
    amino_group_matches = mol.GetSubstructMatches(amino_group_pattern)
    if not amino_group_matches:
        return False, "No primary amino group found"
    else:
        num_amino_groups = len(amino_group_matches)
    
    # Identify peptide bonds
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if peptide_bond_matches:
        return False, "Contains peptide bonds, likely a peptide"
    
    # Check if amino groups are connected to carbon atoms (exclude those attached to heteroatoms)
    free_amino_groups = []
    for match in amino_group_matches:
        nitrogen_idx = match[0]
        atom = mol.GetAtomWithIdx(nitrogen_idx)
        is_attached_to_carbon = any([
            neighbor.GetAtomicNum() == 6
            for neighbor in atom.GetNeighbors()
        ])
        if is_attached_to_carbon:
            free_amino_groups.append(nitrogen_idx)
    
    if not free_amino_groups:
        return False, "No free amino group attached to carbon"
    else:
        num_free_amino_groups = len(free_amino_groups)
    
    # Check for molecular weight to exclude large peptides
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too high for typical amino acids"
    
    return True, f"Contains {num_carboxylic_acids} carboxylic acid group(s) and {num_free_amino_groups} amino group(s)"