"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:35244 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidines have the structure RC(=NR)NR2 and contain the -C(=NH)-NH2 group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Carboxamidine SMARTS pattern
    carboxamidine_pattern = Chem.MolFromSmarts("C(=[N&D2])[N&H2,H1]")
    
    # Check if the molecule contains the carboxamidine group
    matches = mol.GetSubstructMatches(carboxamidine_pattern)
    
    if matches:
        for match in matches:
            # Check if the matched atoms are not part of a ring
            if not any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
                # Check if the nitrogen atoms are not further substituted
                n1_idx, n2_idx = match[1], match[2]
                n1_atom = mol.GetAtomWithIdx(n1_idx)
                n2_atom = mol.GetAtomWithIdx(n2_idx)
                if n1_atom.GetTotalDegree() == 1 and n2_atom.GetTotalDegree() == 2:
                    return True, "Molecule contains the carboxamidine group RC(=NR)NR2"
    
    return False, "No carboxamidine group found"