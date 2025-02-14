"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: CHEBI:32935 aralkylamine
An alkylamine in which the alkyl group is substituted by an aromatic group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Look for aromatic atoms
    aromatic_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if not aromatic_atoms:
        return False, "No aromatic atoms found"
    
    # Look for primary/secondary amine group(s)
    amine_pattern = Chem.MolFromSmarts("[NH2,NH1]")
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    if not amine_matches:
        return False, "No primary or secondary amine groups found"
    
    # Check if amine group is connected to an aromatic system (possibly via an alkyl chain)
    for amine_idx in amine_matches:
        amine_atom = mol.GetAtomWithIdx(amine_idx)
        for aromatic_idx in aromatic_atoms:
            if mol.AreAtomsSeparatedByAnyRing(amine_idx, aromatic_idx):
                continue  # Amine and aromatic atom are on opposite sides of a ring
            path = Chem.rdmolops.GetShortestPath(mol, amine_idx, aromatic_idx)
            if path is not None:
                alkyl_chain = all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path[1:-1])
                if alkyl_chain:
                    return True, "Contains an alkylamine with aromatic substituent"
    
    # Amine group(s) not connected to aromatic system via alkyl chain
    return False, "Amine group(s) not connected to aromatic system via alkyl chain"