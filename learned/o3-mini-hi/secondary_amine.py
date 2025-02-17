"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine 
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
This function attempts to detect a secondary amine group. It checks each nitrogen atom to see if it has two 
hydrocarbyl (carbon) substituents and one hydrogen. In case the nitrogen is nitrosated (i.e. having an N=O substituent)
we treat that nitroso as “virtually restoring” the hydrogen so that N-nitrosamines (e.g. N-nitrosopiperidine) are classified.
We also exclude cases where one of the groups on nitrogen is a carbonyl (i.e. an amide bond) so that amide nitrogens are not mis‐classified.
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains a secondary amine functional group based on its SMILES string.
    A secondary amine should be derived from ammonia by replacing two hydrogens with hydrocarbyl groups.
    In our algorithm we look for a nitrogen atom (atomic number 7) that is attached to:
      • exactly two carbon neighbors (hydrocarbyl substituents) (and ensuring these carbons are not part of a carbonyl group)
      • and an effective one hydrogen.  Here “effective” means that if one neighbor is a nitroso moiety (-N=O), 
        it is treated as replacing a hydrogen.
        
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if a secondary amine is found, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue

        # Retrieve the total number of hydrogens (implicit + explicit)
        h_count = atom.GetTotalNumHs()

        # Counters:
        carbon_count = 0   # count of substituents that are carbon AND are not in a carbonyl bond with this N
        nitroso_count = 0  # count of substituents that are nitroso groups (i.e. an N that is doubly bonded to an O)
        other_count = 0    # any other substituents (such as atoms that should not be counted as hydrocarbyl)

        # Iterate over neighbors of the candidate nitrogen.
        for nbr in atom.GetNeighbors():
            # Get the bond between the nitrogen and this neighbor.
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # Function to decide if a carbon neighbor is in a carbonyl group (C=O) with respect to this bond.
            def is_carbonyl(nbr_atom):
                # Only consider carbon atoms.
                if nbr_atom.GetAtomicNum() != 6:
                    return False
                # Check all bonds on this carbon (except the one to our nitrogen) to see if there is any double bond to oxygen.
                for bond2 in nbr_atom.GetBonds():
                    # Skip the bond back to our nitrogen.
                    if bond2.GetOtherAtom(nbr_atom).GetIdx() == atom.GetIdx():
                        continue
                    # If a double bond to O is found, treat this substituent as a carbonyl.
                    if bond2.GetBondTypeAsDouble() == 2 and bond2.GetOtherAtom(nbr_atom).GetAtomicNum() == 8:
                        return True
                return False

            # If neighbor is a carbon and not carbonyl then consider it a hydrocarbyl group.
            if nbr.GetAtomicNum() == 6:
                if is_carbonyl(nbr):
                    other_count += 1  # flag as not a simple hydrocarbyl substituent.
                else:
                    carbon_count += 1
            # If neighbor is a nitrogen, check if it is part of a nitroso group (i.e. has a double bond to oxygen)
            elif nbr.GetAtomicNum() == 7:
                # Look for a double bond from this neighbor to an oxygen.
                has_double_o = False
                for subnbr in nbr.GetNeighbors():
                    # Skip if subnbr is our original atom (avoid backtracking).
                    if subnbr.GetIdx() == atom.GetIdx():
                        continue
                    if subnbr.GetAtomicNum() == 8:
                        bond2 = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        if bond2 and bond2.GetBondTypeAsDouble() == 2:
                            has_double_o = True
                            break
                if has_double_o:
                    nitroso_count += 1
                else:
                    other_count += 1
            else:
                # Any substituent that is not carbon or (nitroso) nitrogen is flagged as 'other'
                other_count += 1

        # Define effective hydrogen count: if a nitroso group is present, treat it as replacing a hydrogen.
        effective_h = h_count + nitroso_count

        # We now check that the candidate nitrogen has exactly two hydrocarbyl (carbon) substituents, no disallowed other groups,
        # and an effective hydrogen count of 1.
        if carbon_count == 2 and other_count == 0 and effective_h == 1:
            return True, "Found a nitrogen with exactly two carbon substituents and one effective hydrogen (including nitroso correction)"
    
    return False, "No secondary amine found: no nitrogen atom with exactly two appropriate hydrocarbyl substituents and one hydrogen"

# Example usage (for debugging):
if __name__ == "__main__":
    test_smiles = [
        "CNc1ccccc1",               # N-methylaniline, a secondary amine.
        "O=NN1CCCCC1",              # N-nitrosopiperidine: should be classified as secondary amine.
        "[H]N(C)C",                # dimethylamine, a secondary amine.
        "O=C(NCc1ccccc1)C"          # an amide: should not be classified.
    ]
    for smi in test_smiles:
        result, reason = is_secondary_amine(smi)
        print(f"SMILES: {smi} -> {result}, Reason: {reason}")