"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N-acylglycine
An N-acyl-amino acid in which amino acid specified is glycine.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    N-acylglycines have a glycine moiety (NH-CH2-COOH) with an acyl group (R-C=O-) 
    attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Convert to neutral form if possible (handle cases with charged carboxylate)
    mol = Chem.AddHs(mol)

    # Look for glycine moiety pattern (NH-CH2-COOH)
    # [NH;!$(NC=O):1] means nitrogen with H, not already part of an amide
    # [CH2:2] means CH2 group
    # [CX3:3](=[OX1])[OH] means carboxylic acid group
    glycine_pattern = Chem.MolFromSmarts('[NH;!$(NC=O):1][CH2:2][CX3:3](=[OX1])[OH]')
    
    # Look for acyl group pattern (R-C(=O)-)
    acyl_pattern = Chem.MolFromSmarts('[#6,#1;A][CX3](=O)[NX3;H1]')
    
    glycine_matches = mol.GetSubstructMatches(glycine_pattern)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    
    if not glycine_matches:
        return False, "No glycine moiety (NH-CH2-COOH) found"
        
    if not acyl_matches:
        return False, "No acyl group (R-C=O-NH) found"
    
    # Now verify that the acyl group is connected to the glycine nitrogen
    # Get all atoms in the molecule
    atoms = mol.GetAtoms()
    
    # Check each glycine match
    for gly_match in glycine_matches:
        gly_n = mol.GetAtomWithIdx(gly_match[0])  # nitrogen atom
        gly_c = mol.GetAtomWithIdx(gly_match[1])  # CH2 carbon
        gly_acid = mol.GetAtomWithIdx(gly_match[2])  # acid carbon
        
        # Check neighbors of nitrogen
        for neighbor in gly_n.GetNeighbors():
            # Look for carbonyl carbon attached to nitrogen
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for n_neighbor in neighbor.GetNeighbors():
                    if n_neighbor.GetAtomicNum() == 8 and n_neighbor.GetBonds()[0].GetBondType() == Chem.BondType.DOUBLE:
                        # Found N-acylglycine structure
                        return True, "Contains glycine moiety with N-acyl group"
    
    return False, "Glycine and acyl groups present but not connected correctly"