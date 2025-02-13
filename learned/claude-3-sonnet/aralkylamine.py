"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: CHEBI:35484 aralkylamine

An aralkylamine is an alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for disconnected structures
    mol_frags = Chem.GetMolFrags(mol, aromatics=Chem.SMARTSFromSmarts('c'))
    
    for frag in mol_frags:
        frag_mol = Chem.MolFragmentToSmiles(mol, frag)
        
        # Look for alkylamino group (-N-) and aromatic ring
        alkylamine_pattern = Chem.MolFromSmarts("[N;!$(NC=O)!$(N([#6])[#6])]")
        aromatic_pattern = Chem.MolFromSmarts("a")
        
        frag_mol = Chem.MolFromSmiles(frag_mol)
        alkylamine_matches = frag_mol.GetSubstructMatches(alkylamine_pattern)
        aromatic_matches = frag_mol.GetSubstructMatches(aromatic_pattern)
        
        if alkylamine_matches and aromatic_matches:
            # Check for reasonable molecular weight and atom count
            mol_wt = rdMolDescriptors.CalcExactMolWt(frag_mol)
            if mol_wt > 100 and mol_wt < 500:
                n_atoms = frag_mol.GetNumAtoms()
                if n_atoms > 10 and n_atoms < 50:
                    return True, "Contains alkylamino group and aromatic ring"
    
    return False, "No alkylamino group or aromatic ring found"