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
    
    # Look for alkylamino group (-N-) and aromatic ring attached to alkyl chain
    alkylamine_pattern = Chem.MolFromSmarts("[N;!$(NC=O)]")
    aromatic_pattern = Chem.MolFromSmarts("a")
    
    alkylamine_matches = mol.GetSubstructMatches(alkylamine_pattern)
    aromatic_matches = mol.GetSubstructMatches(aromatic_pattern)
    
    if not alkylamine_matches or not aromatic_matches:
        return False, "No alkylamino group or aromatic ring found"
    
    # Check if aromatic ring is directly attached to alkyl chain of alkylamino group
    for amine_idx in alkylamine_matches:
        amine_atom = mol.GetAtomWithIdx(amine_idx)
        for neighbor in amine_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                # Check for alkyl chain between amine and aromatic ring
                alkyl_chain = Chem.Mol(Chem.FindEqualInsatNEndPts(Chem.MolFromSmiles(Chem.MolToSmiles(Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(mol, amine_idx, 4), atomIdxList=[amine_idx])), allNodes=True))
                if alkyl_chain.GetNumAtoms() > 1:
                    return True, "Contains alkylamino group with aromatic ring attached via alkyl chain"
    
    return False, "Aromatic ring not directly attached to alkyl chain of alkylamino group"