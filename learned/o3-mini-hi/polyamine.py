"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.
"""

from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    
    The algorithm works as follows:
      1. The SMILES is parsed into an RDKit molecule.
      2. For each nitrogen atom (atomic number 7) in the molecule, we check:
           - It is not involved in an amide group (i.e. directly bonded to a carbon that is part of a carbonyl, C(=O)).
           - If the nitrogen is part of a ring, we require that it have at least one hydrogen (to avoid counting embedded heterocyclic nitrogens).
           - It is not quaternary (i.e. has fewer than four neighbors).
      3. If at least two such amino groups are found, then the molecule is classified as a polyamine.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polyamine, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    amino_count = 0  # counter for amino groups
    
    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # skip non-nitrogen atoms
        
        # Check if this nitrogen is part of an amide.
        # We define a nitrogen as “amide-like” if it is directly bonded to a carbon that
        # is double-bonded to an oxygen.
        is_amide = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # check for carbon neighbor
                for bond in neighbor.GetBonds():
                    # Skip the bond that connects neighbor and our nitrogen
                    if bond.GetOtherAtom(neighbor).GetIdx() == atom.GetIdx():
                        continue
                    # If neighbor carbon is double-bonded to an oxygen, mark as amide.
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(neighbor).GetAtomicNum() == 8:
                        is_amide = True
                        break
                if is_amide:
                    break
        if is_amide:
            continue  # skip this nitrogen; it is part of an amide
        
        # For ring nitrogen atoms, we only consider those bearing at least one hydrogen.
        # (This helps avoid counting heterocyclic nitrogen atoms that do not act as free amino groups)
        if atom.IsInRing() and atom.GetTotalNumHs() == 0:
            continue
        
        # Exclude quaternary ammonium centers (4 substituents)
        if atom.GetDegree() >= 4:
            continue
        
        # If we reached here, we consider this nitrogen as an amino group.
        amino_count += 1

    # Check if we found at least two amino groups.
    if amino_count >= 2:
        return True, f"Contains {amino_count} amino groups."
    else:
        return False, f"Found only {amino_count} amino group(s), need at least 2 to be a polyamine."
        
# Example usage:
if __name__ == "__main__":
    # Examples provided (you can try one or more)
    test_smiles = [
        "CCNc1nc(N)nc(O)n1",           # 4-amino-6-(ethylamino)-1,3,5-triazin-2-ol
        "NCCCNCCNCCCN",                # 3,2,3-tetramine
        "NCCN",                        # ethylenediamine
        "NCN"                          # methanediamine
    ]
    for s in test_smiles:
        result, reason = is_polyamine(s)
        print(f"SMILES: {s} -> {result} ({reason})")