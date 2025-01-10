"""
Classifies: CHEBI:27096 trichlorobenzene
"""
from rdkit import Chem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule contains a trichlorobenzene motif based on its SMILES string.
    A trichlorobenzene motif is defined as a benzene ring with exactly three chlorine atom substituents,
    potentially part of larger molecules.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains a trichlorobenzene motif, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all aromatic rings in the molecule
    aromatic_rings = mol.GetAromaticRings()
    
    # For each aromatic ring, count the number of chlorine substituents
    for ring in aromatic_rings:
        chlorine_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 17:  # Check if atom is Chlorine (Cl)
                chlorine_count += 1
        
        # Check if the ring contains exactly 3 chlorine atoms
        if chlorine_count == 3:
            return True, "Contains a trichlorobenzene motif (benzene ring with three chlorine atoms)"

    return False, "Does not contain a trichlorobenzene motif"
    
# Example usage
smiles_strings = [
    "Clc1cc(Cl)cc(Cl)c1",     # 1,3,5-trichlorobenzene
    "Clc1cccc(Cl)c1Cl",       # 1,2,3-trichlorobenzene
    "Clc1ccc(Cl)c(Cl)c1",     # 1,2,4-trichlorobenzene
    "COc1cc(Cl)c(Cl)cc1Cl"    # 2,4,6-trichloroanisole
]

for smiles in smiles_strings:
    result, reason = is_trichlorobenzene(smiles)
    print(f"SMILES: {smiles} => Result: {result}, Reason: {reason}")