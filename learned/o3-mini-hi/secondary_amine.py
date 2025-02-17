"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine 
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
This improved function looks for a nitrogen (atomic number 7) that is bonded to exactly two carbon (hydrocarbyl) groups and one effective hydrogen.
A nitroso substituent (-N=O) is treated as replacing a hydrogen. Additionally, we exclude candidate nitrogens that are directly bonded to a
carbonyl group (i.e. are part of an amide), because peptide and similar acylated nitrogens, although they might match the simple count,
should be excluded.
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains a secondary amine functional group based on its SMILES string.
    A secondary amine should be derived from ammonia by replacing two hydrogens with hydrocarbyl groups.
    
    The algorithm:
      1. For each nitrogen atom (atomic number 7) in the molecule, first reject it if it is part of an amide linkage.
         We do this by checking if any carbon neighbor is a carbonyl carbon (bonded by a double bond to an oxygen).
      2. For each remaining candidate N, count:
           - the total number of hydrogens (implicit+explicit)
           - the number of carbon (atomic number 6) substituents that are NOT carbonyl carbons
           - the number of nitroso substituents, defined as a neighboring nitrogen that is doubly bonded to an oxygen.
      3. Define the effective hydrogen count as the sum of the hydrogen count and the nitroso count.
         If a candidate N has exactly 2 carbon substituents (non-carbonyl) and an effective hydrogen count of 1, we flag it as a secondary amine.
         
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if a secondary amine is found, False otherwise.
       str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    def is_amide_nitrogen(atom, mol):
        """
        Checks if the nitrogen atom is part of an amide linkage.
        That is, if any carbon neighbor is directly bonded to an oxygen by a double bond.
        """
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                # Obtain the bond between the nitrogen and its carbon neighbor.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # Check all bonds of the carbon neighbor (skipping the bond to the nitrogen)
                for bond2 in nbr.GetBonds():
                    other = bond2.GetOtherAtom(nbr)
                    if other.GetIdx() == atom.GetIdx():
                        continue
                    # If there is a double bond to oxygen, the carbon is carbonyl.
                    if other.GetAtomicNum() == 8 and bond2.GetBondType() == Chem.BondType.DOUBLE:
                        return True
        return False

    # Loop over every nitrogen atom in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue

        # Exclude if nitrogen is in an amide
        if is_amide_nitrogen(atom, mol):
            continue

        # Count total hydrogens (implicit + explicit)
        h_count = atom.GetTotalNumHs()
        carbon_count = 0   # Count of carbon substituents that are simple hydrocarbyl (not part of a carbonyl)
        nitroso_count = 0  # Count of substituents that are nitroso groups.
        other_count = 0    # Any other substituents.

        # Iterate over neighbors of the candidate nitrogen.
        for nbr in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # For carbon atoms, check if the bond is to a carbonyl.
            if nbr.GetAtomicNum() == 6:
                # Check if this carbon acts as a carbonyl partner.
                is_carbonyl = False
                for bond2 in nbr.GetBonds():
                    other = bond2.GetOtherAtom(nbr)
                    # Skip the bond back to the nitrogen.
                    if other.GetIdx() == atom.GetIdx():
                        continue
                    if other.GetAtomicNum() == 8 and bond2.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl = True
                        break
                if is_carbonyl:
                    # Flag as an undesired substituent since it implies an acyl group.
                    other_count += 1
                else:
                    carbon_count += 1

            # For nitrogen neighbors, check if they are part of a nitroso group.
            elif nbr.GetAtomicNum() == 7:
                has_double_o = False
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetIdx() == atom.GetIdx():
                        continue
                    if subnbr.GetAtomicNum() == 8:
                        bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bond2 and bond2.GetBondType() == Chem.BondType.DOUBLE:
                            has_double_o = True
                            break
                if has_double_o:
                    nitroso_count += 1
                else:
                    other_count += 1
            else:
                other_count += 1

        # Effective hydrogen count: count the hydrogens plus any nitroso substituent (which is treated as replacing a hydrogen).
        effective_h = h_count + nitroso_count

        # Check for secondary amine: exactly two non-carbonyl carbon substituents, no extra groups,
        # and an effective hydrogen count of one.
        if carbon_count == 2 and other_count == 0 and effective_h == 1:
            return True, ("Found a nitrogen with exactly two carbon substituents and one effective hydrogen "
                          "(including nitroso correction) and no amide character")
    
    return False, "No secondary amine found: no nitrogen atom with exactly two appropriate hydrocarbyl substituents and one hydrogen"

# Example usage (for debugging):
if __name__ == "__main__":
    test_smiles = [
        "C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O",  # (R)-dobutamine: expects True
        "CNC1CCCCC1",  # N-methylcyclohexylamine: True
        "CC(C)NCC(C)(C)N",  # N(1)-isopropyl-2-methylpropan-1,2-diamine: True
        "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",  # (S)-dobutamine: True
        "[H]C(=O)CCCNCCCN",  # N-(3-aminopropyl)-4-aminobutanal: True
        "CNc1ccccc1",  # N-methylaniline: True
        "O=NN1CCCCC1",  # N-nitrosopiperidine: True
        "[H]N(C)C",  # dimethylamine: True
        "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1",  # A peptide fragment: should be False (amide nitrogen)
    ]
    for smi in test_smiles:
        result, reason = is_secondary_amine(smi)
        print(f"SMILES: {smi}\n   -> {result}, Reason: {reason}\n")