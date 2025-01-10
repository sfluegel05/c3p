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

    # Check for carboxyl group (-COOH) attached to the alpha-carbon
    carboxyl_pattern = Chem.MolFromSmarts("[CX4H]([CX3](=[OX1])[OX2H1])")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group attached to alpha-carbon"

    # Check for amino group (-NH2) attached to the alpha-carbon
    amino_pattern = Chem.MolFromSmarts("[CX4H]([NX3H2])")
    if not mol.HasSubstructMatch(amino_pattern):
        return False, "No amino group attached to alpha-carbon"

    # Check for L-configuration at the alpha-carbon
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return False, "No chiral center found"
    
    # Verify that the alpha-carbon has the L-configuration
    alpha_carbon = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 4:  # Potential alpha-carbon
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'N' in neighbors and 'C' in neighbors:  # Amino and carboxyl groups
                alpha_carbon = atom.GetIdx()
                break
    
    if alpha_carbon is None:
        return False, "No alpha-carbon found"
    
    # Check if the alpha-carbon is chiral and has L-configuration
    for center, config in chiral_centers:
        if center == alpha_carbon:
            if config == 'R':  # L-configuration corresponds to S in RDKit
                return False, "Alpha-carbon has D-configuration"
            elif config == 'S':
                return True, "Alpha-carbon has L-configuration"
            else:
                return False, "Chirality at alpha-carbon is unassigned"
    
    return False, "No L-configuration found at alpha-carbon"