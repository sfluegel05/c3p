"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    
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
    
    # Identify all chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    # Attempt to find alpha-carbon with L-configuration
    for center, chirality in chiral_centers:
        atom = mol.GetAtomWithIdx(center)
        if atom.GetAtomicNum() == 6:  # Carbon
            # Make sure we have an NH group and a carboxylic acid (COOH)
            n_neighbor = False
            c_group = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 7 and nbr.GetTotalNumHs() > 0:
                    n_neighbor = True
                if nbr.GetAtomicNum() == 8:
                    o_neighbors = [atom.GetAtomicNum() for atom in nbr.GetNeighbors()]
                    if o_neighbors.count(8) == 1:  # Connected to another oxygen
                        c_group = True
            
            # Check if chiral center has 'S' configuration (commonly corresponds to L-configuration)
            if n_neighbor and c_group and chirality == 'S':
                return True, "Chiral center with L-configuration, amino and carboxyl groups identified"
    
    return False, "No chiral center suitable for L-alpha-amino acid found"