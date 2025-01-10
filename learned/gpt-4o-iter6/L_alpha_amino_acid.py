"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Identify chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    
    # Define SMARTS pattern for an L-alpha-amino acid backbone (N[C@@H](C)C(O)=O)
    aa_smarts = '[C@H](N)C(=O)O'
    aa_pattern = Chem.MolFromSmarts(aa_smarts)
    
    if not mol.HasSubstructMatch(aa_pattern):
        return False, "No L-alpha-amino acid backbone pattern found"
    
    # Check if there's at least one 'S' configuration chiral center (common for L-amino acids)
    for center, chirality in chiral_centers:
        if chirality in ['S', '(S)']:
            atom = mol.GetAtomWithIdx(center)
            if atom.GetSymbol() == 'C':
                # Verify the core structure and known neighbours
                n_neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
                o_count = n_neighbors.count('O')
                n_count = n_neighbors.count('N')
                
                if n_count == 1 and o_count >= 1:
                    return True, "Chiral center with L-configuration identified, amino and carboxyl groups confirmed."
    
    return False, "No suitable chiral configuration or functional groups for L-alpha-amino acid found."