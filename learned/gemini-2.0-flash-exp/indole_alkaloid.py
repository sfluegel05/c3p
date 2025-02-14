"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid contains an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More robust SMARTS pattern for indole - simplified to a 5-membered ring with N connected to an aromatic 6-membered ring.
    indole_pattern = Chem.MolFromSmarts("c1cc2[nH]c(c1)cc2") # simplified smarts, aromatic C connected to aromatic N

    if mol.HasSubstructMatch(indole_pattern):
        
        #Verify the substructure has at least 9 ring atoms and 1 Nitrogen.
        match = mol.GetSubstructMatch(indole_pattern)
        ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
        
        n_count = 0
        aromatic_count = 0
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 7: #nitrogen
                n_count+=1
            if atom.GetIsAromatic():
              aromatic_count +=1
        
        if n_count >0 and aromatic_count>=9:
            return True, "Contains an indole substructure"

    return False, "Does not contain an indole substructure"