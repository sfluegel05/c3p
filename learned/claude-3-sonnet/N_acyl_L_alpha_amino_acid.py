"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible pattern for alpha-amino acid core
    # Matches both protonated and deprotonated forms
    # Uses [C@@H] or [C@H] for L-amino acids
    amino_acid_pattern = Chem.MolFromSmarts('[C@@H,C@H]([NH0,NH1,NH2])(C(=O)[OH0,OH1,O-])[#6]')
    
    if not mol.HasSubstructMatch(amino_acid_pattern, useChirality=True):
        return False, "No L-alpha-amino acid core found"

    # More flexible N-acyl pattern
    # Matches any acyl group (R-C=O) attached to the amino acid nitrogen
    n_acyl_pattern = Chem.MolFromSmarts('[#6]C(=O)[NH0,NH1][C@@H,C@H]([#6])[CX3](=[OX1])[OH0,OH1,O-]')
    
    if not mol.HasSubstructMatch(n_acyl_pattern, useChirality=True):
        return False, "No N-acyl group found"

    # Get matches of the amino acid pattern
    matches = mol.GetSubstructMatches(amino_acid_pattern, useChirality=True)
    
    for match in matches:
        alpha_carbon_idx = match[0]
        nitrogen_idx = match[1]
        
        # Get the alpha carbon
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        
        # Skip if chirality is not specified
        if alpha_carbon.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
            
        # Get the nitrogen atom
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        
        # Check nitrogen's acylation
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                # Look for carbonyl group attached to nitrogen
                carbonyl_found = False
                for n2 in neighbor.GetNeighbors():
                    if n2.GetAtomicNum() == 8:  # Oxygen
                        if any(b.GetBondType() == Chem.BondType.DOUBLE 
                              for b in neighbor.GetBonds()):
                            carbonyl_found = True
                            break
                
                if carbonyl_found:
                    # Verify carboxylic acid/carboxylate group
                    for c_neighbor in alpha_carbon.GetNeighbors():
                        if c_neighbor.GetAtomicNum() == 6:  # Carbon
                            if any(n.GetAtomicNum() == 8 for n in c_neighbor.GetNeighbors()):
                                if any(b.GetBondType() == Chem.BondType.DOUBLE 
                                      for b in c_neighbor.GetBonds()):
                                    return True, "Contains N-acyl-L-alpha-amino acid structure"

    return False, "No valid N-acyl-L-alpha-amino acid structure found"