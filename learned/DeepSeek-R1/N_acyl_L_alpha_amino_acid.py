"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    Criteria: L-alpha-amino acid structure with at least one N-acyl (amide) substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for at least one N-acyl group (amide)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No N-acyl group found"

    # Find potential alpha carbons (connected to carboxyl group)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Check each carboxyl group's alpha carbon for amino acid structure
    for match in carboxyl_matches:
        carboxyl_c = match[0]
        alpha_c = None
        
        # Find alpha carbon (adjacent to carboxyl carbon)
        for neighbor in mol.GetAtomWithIdx(carboxyl_c).GetNeighbors():
            if neighbor.GetSymbol() == 'C':
                alpha_c = neighbor.GetIdx()
                break
        if alpha_c is None:
            continue

        # Verify alpha carbon has nitrogen connection (amino/amide)
        has_nitrogen = any(atom.GetSymbol() == 'N' 
                          for atom in mol.GetAtomWithIdx(alpha_c).GetNeighbors())
        if not has_nitrogen:
            continue

        # Check stereochemistry (L-configuration = S)
        mol_h = Chem.AddHs(mol)
        chiral_centers = Chem.FindMolChiralCenters(mol_h)
        for center, chirality in chiral_centers:
            if center == alpha_c and chirality == 'S':
                return True, "L-alpha-amino acid with N-acyl substituent"

    return False, "Not an L-alpha-amino acid or missing N-acyl group"