"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:35244 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidines have the structure RC(=NR)NR2 and contain the -C(=N)-NH2 group.

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
    
    # Generate tautomers and resonance structures
    tautomers = list(AllChem.EnumerateTautomers(mol))
    
    # Carboxamidine SMARTS pattern
    carboxamidine_pattern = Chem.MolFromSmarts("C(=[N&D2])N")
    
    # Check if any tautomer contains the carboxamidine group
    for tautomer in tautomers:
        if tautomer.HasSubstructMatch(carboxamidine_pattern):
            # Check if the matched group is not part of a ring
            matched_atoms = tautomer.GetSubstructMatches(carboxamidine_pattern)[0]
            if not tautomer.GetAtomWithIdx(matched_atoms[0]).IsInRing():
                return True, "Molecule contains the carboxamidine group RC(=NR)NR2"
    
    return False, "No carboxamidine group found"