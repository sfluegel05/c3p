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
    
    # Identify all carbon atoms with 4 single bonds (sp3 hybridized)
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    
    alpha_carbon = None
    for center, chirality in chiral_centers:
        atom = mol.GetAtomWithIdx(center)
        # Ensure the carbon has one -NH2, one -COOH, and possibly other groups
        if atom.GetAtomicNum() == 6:  # Carbon
            # Check for amino (NH2) and carboxylic acid (COOH) attachments
            neighbors = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
            if (neighbors.count(7) == 1 and  # One nitrogen neighbor
                neighbors.count(8) >= 1):    # At least one oxygen neighbor
                # Verify this is the alpha carbon
                n_bound, o_bound = False, False
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 7:
                        # Check if nitrogen has 2 hydrogens (NH2)
                        n_bound = sorted(nbr.GetNeighbors(), key=lambda x: x.GetAtomicNum())[-1].GetAtomicNum() == 1
                    if nbr.GetAtomicNum() == 8:
                        # Check if there are two oxygens in carboxylate formation (COO)
                        num_o = sum(n.GetAtomicNum() == 8 for n in nbr.GetNeighbors())
                        if num_o == 2:
                            o_bound = True
                
                if n_bound and o_bound:
                    if chirality == 'S':
                        alpha_carbon = center
                        break
    
    if alpha_carbon is None:
        return False, "No chiral center suitable for L-alpha-amino acid found"

    return True, "Chiral center with L-configuration, amino and carboxyl groups identified"