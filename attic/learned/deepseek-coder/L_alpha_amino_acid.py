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

    # More robust carboxyl group pattern (including zwitterionic forms)
    carboxyl_patterns = [
        Chem.MolFromSmarts("[CX3](=O)[OX1H0-,OX2H1]"),  # Normal carboxyl
        Chem.MolFromSmarts("[CX3](=O)[O-]"),           # Zwitterionic form
        Chem.MolFromSmarts("[CX3](=O)O")               # Protonated form
    ]
    
    carboxyl_matches = []
    for pattern in carboxyl_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            carboxyl_matches.extend(matches)
    
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # More robust amino group pattern (including zwitterionic forms)
    amino_patterns = [
        Chem.MolFromSmarts("[NX3;H2,H1,H0]"),          # Normal amino
        Chem.MolFromSmarts("[NH3+]"),                  # Zwitterionic form
        Chem.MolFromSmarts("[NX3H2]")                  # Protonated form
    ]
    
    amino_matches = []
    for pattern in amino_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            amino_matches.extend(matches)
    
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

    alpha_carbon = None
    for amino_match in amino_matches:
        amino_nitrogen = amino_match[0]
        amino_neighbors = mol.GetAtomWithIdx(amino_nitrogen).GetNeighbors()
        for neighbor in amino_neighbors:
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() in alpha_candidates:
                alpha_carbon = neighbor.GetIdx()
                break
        if alpha_carbon is not None:
            break

    if alpha_carbon is None:
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