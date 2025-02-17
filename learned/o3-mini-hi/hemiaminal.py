"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemi-aminal compounds
Definition: Any organic amino compound that has an amino group and a hydroxy group 
attached to the same tetrahedral (sp3) carbon that carries exactly one hydrogen.
Hemiaminals are intermediates in imine formation.
Improvement: Instead of strictly requiring exactly three heavy atom bonds, we 
simply require that an sp3 (non-aromatic) carbon with exactly one hydrogen has at 
least one single-bonded OH group (oxygen with ≥1 hydrogen) and at least one single-
bonded amino group (nitrogen). This "looser" matching may catch hemiaminals embedded
in complex architectures while hopefully reducing false negatives.
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule contains a hemiaminal motif.
    A hemiaminal is defined here as a tetrahedral (sp3, non-aromatic) carbon 
    carrying exactly one hydrogen that is connected via single bonds to both a
    hydroxyl group (-OH) (where the oxygen itself has at least one hydrogen) and an 
    amino group (-NH, -NHR, or -NR2).

    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule contains a hemiaminal motif, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms to search for candidate hemiaminal carbons.
    for atom in mol.GetAtoms():
        # Focus on carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        # Require sp3 hybridization and non-aromatic carbon.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or atom.GetIsAromatic():
            continue
        
        # We require that the carbon carries exactly one hydrogen (implicit plus explicit).
        if atom.GetTotalNumHs() != 1:
            continue

        # Flags for hydroxyl and amino substituents.
        found_OH = False
        found_NH = False

        # Check each neighbor; allow additional substituents.
        for neighbor in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            # Only consider single bonds.
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue

            # Check for hydroxyl group: neighbor is oxygen and has at least one H.
            if neighbor.GetAtomicNum() == 8:
                # For an -OH group, ensure the oxygen has at least one hydrogen, 
                # either explicit or implicit.
                if neighbor.GetTotalNumHs() >= 1:
                    found_OH = True

            # Check for amine group: neighbor is nitrogen.
            elif neighbor.GetAtomicNum() == 7:
                found_NH = True

        # If both -OH and -NH groups are attached, we consider the motif found.
        if found_OH and found_NH:
            reason = ("Molecule contains a hemiaminal motif: a tetrahedral, non-aromatic sp3 carbon with exactly one hydrogen "
                      "that is directly bound (via single bonds) to both a hydroxyl group (-OH) and an amino group (-NH, -NHR, or -NR2).")
            return True, reason

    return False, "No hemiaminal motif (an sp3 carbon with one hydrogen bound to both -OH and -NH groups) was found."

# Example usage (for testing, uncomment the following lines):
# test_smiles_list = [
#     "NC(O)C(O)=O",  # alpha-hydroxyglycine, should be positive.
#     "C1CCCCC1"     # Cyclohexane, should be negative.
# ]
# for smi in test_smiles_list:
#     result, explanation = is_hemiaminal(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")