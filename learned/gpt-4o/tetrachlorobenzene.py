"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a benzene ring carrying four chlorine groups at unspecified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a benzene ring with any four chlorines
    tetrachloro_pattern = Chem.MolFromSmarts("c1(c(cccc1)Cl)Cl.Cl.Cl.Cl")
    
    # Check for substructure match with SMARTS pattern
    if mol.HasSubstructMatch(tetrachloro_pattern):
        # Verify the existence of a benzene ring fulfilling the criteria
        for ring in mol.GetRingInfo().AtomRings():
            if len(ring) == 6:  # Ensure it's a benzene ring
                chloro_count = sum(1 for atom_idx in ring 
                                   if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 17)
                if chloro_count == 4:
                    return True, "Contains a benzene ring with four chlorine atoms"
    
    return False, "Does not match the tetrachlorobenzene pattern"