"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:33762 secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a carbonyl group and a hydroxy group linked
    by a carbon bearing one hydrogen and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and generate resonance structures
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.RemoveHs(mol)
    res_mols = list(AllChem.ResonanceEnumerator(mol))
    
    # Define SMARTS patterns
    carbonyl_smarts = "[C=O]"
    hydroxy_smarts = "[OX2H]"
    alpha_carbon_smarts = "[CX3H1]=[OX1]"
    
    for res_mol in res_mols:
        # Find carbonyl and hydroxy groups
        carbonyl_atoms = [atom.GetIdx() for atom in res_mol.GetAtoms() if atom.HasQuery(Chem.MolFromSmarts(carbonyl_smarts))]
        hydroxy_atoms = [atom.GetIdx() for atom in res_mol.GetAtoms() if atom.HasQuery(Chem.MolFromSmarts(hydroxy_smarts))]
        
        # Check for at least one carbonyl and hydroxy group
        if not carbonyl_atoms or not hydroxy_atoms:
            continue
        
        # Find alpha-carbon atoms
        alpha_carbon_atoms = [atom.GetIdx() for atom in res_mol.GetAtoms() if atom.HasQuery(Chem.MolFromSmarts(alpha_carbon_smarts))]
        
        # Check if any alpha-carbon atom is linked to both carbonyl and hydroxy group
        for alpha_idx in alpha_carbon_atoms:
            alpha_atom = res_mol.GetAtomWithIdx(alpha_idx)
            neighbors = [res_mol.GetAtomWithIdx(n.GetIdx()) for n in alpha_atom.GetNeighbors()]
            if any(n.HasQuery(Chem.MolFromSmarts(carbonyl_smarts)) for n in neighbors) and any(n.HasQuery(Chem.MolFromSmarts(hydroxy_smarts)) for n in neighbors):
                return True, "Contains a secondary alpha-hydroxy ketone moiety"
    
    return False, "No secondary alpha-hydroxy ketone moiety found"