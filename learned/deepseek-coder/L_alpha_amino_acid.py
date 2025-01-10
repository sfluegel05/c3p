"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid has a carboxyl group, an amino group, and L-configuration at the alpha-carbon.

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

    # More flexible carboxyl group pattern (-COO or -COOH)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0,O-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # More flexible amino group pattern (NH2, NH, or N with any substituents)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group found"

    # Find potential alpha-carbons (carbons connected to both amino and carboxyl groups)
    alpha_candidates = set()
    for carboxyl_match in carboxyl_matches:
        carboxyl_carbon = carboxyl_match[0]
        carboxyl_neighbors = mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors()
        for neighbor in carboxyl_neighbors:
            if neighbor.GetSymbol() == 'C':
                alpha_candidates.add(neighbor.GetIdx())

    for amino_match in amino_matches:
        amino_nitrogen = amino_match[0]
        amino_neighbors = mol.GetAtomWithIdx(amino_nitrogen).GetNeighbors()
        for neighbor in amino_neighbors:
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() in alpha_candidates:
                alpha_carbon = neighbor.GetIdx()
                break
        else:
            continue
        break
    else:
        return False, "No alpha-carbon found connecting amino and carboxyl groups"

    # Check chirality at alpha-carbon
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral centers found"

    # Verify L-configuration at alpha-carbon
    for center, config in chiral_centers:
        if center == alpha_carbon:
            if config == 'R':  # L-configuration corresponds to S in RDKit
                return False, "Alpha-carbon has D-configuration"
            elif config == 'S':
                return True, "Alpha-carbon has L-configuration"
            else:
                return False, "Chirality at alpha-carbon is unassigned"
    
    return False, "No L-configuration found at alpha-carbon"