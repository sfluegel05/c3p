"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    L-amino acids have S configuration at the alpha carbon (except cysteine and derivatives).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic alpha-amino acid pattern:
    # [N] - [C@H|@@H] - C(=O)O
    # The pattern looks for a carbon with exactly one hydrogen (alpha carbon)
    # connected to an amine group and a carboxyl group
    aa_pattern = Chem.MolFromSmarts('[NX3;H2,H3,H4+][CX4;H1]([#6,#1])[CX3](=[OX1])[OX2H,OX1-]')
    if not mol.HasSubstructMatch(aa_pattern):
        return False, "No alpha-amino acid pattern found"
    
    matches = mol.GetSubstructMatches(aa_pattern)
    
    # Check each match for correct stereochemistry
    for match in matches:
        alpha_carbon_idx = match[1]  # Second atom in pattern is alpha carbon
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        
        # Must have specified stereochemistry
        if alpha_carbon.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue
            
        # Get canonical SMILES with stereochemistry information
        canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        
        # Find the chirality marker for this carbon in the SMILES
        for i, atom in enumerate(mol.GetAtoms()):
            if atom.GetIdx() == alpha_carbon_idx:
                # Check if this is a cysteine or cysteine derivative
                is_cysteine = False
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Carbon
                        for next_neighbor in neighbor.GetNeighbors():
                            if next_neighbor.GetAtomicNum() == 16:  # Sulfur
                                is_cysteine = True
                                break
                
                # For L-amino acids:
                # - Normal case: S configuration (@@ in SMILES)
                # - Cysteine case: R configuration (@ in SMILES)
                if (is_cysteine and '@' in canonical_smiles and '@@' not in canonical_smiles) or \
                   (not is_cysteine and '@@' in canonical_smiles):
                    return True, "Found L-alpha-amino acid with correct stereochemistry"
                
    # If we get here, we found the pattern but not the correct stereochemistry
    return False, "Found alpha-amino acid pattern but incorrect or unspecified stereochemistry"