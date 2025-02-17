"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated Fatty Acid

Definition: Any fatty acid with one double or triple bond in the fatty acid chain 
(and all remaining carbon–carbon bonds are single). Fatty acids are characterized 
by having a carboxylic acid group, and the monounsaturated ones have exactly one 
non–acid-group carbon–carbon unsaturation. This version tries to restrict the unsaturation 
to occur in the long, acyclic fatty acid chain.
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    
    The function performs these steps:
      1. Parses the SMILES string.
      2. Rejects cyclic molecules (most fatty acid chains are acyclic).
      3. Checks for a carboxylic acid group (recognizing both protonated and deprotonated forms).
      4. Verifies that the acid carbon (of the –COOH group) is terminal (has exactly one carbon neighbor).
      5. Requires that the molecule contains a sufficient number of carbon atoms.
      6. Counts the unsaturation (double or triple C–C bonds) in the molecule
         while ignoring the carbonyl (C=O) bond of the acid group and bonds that are not “internal”
         to the chain (i.e. if one end of the unsaturation is effectively terminal, it is not counted).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monounsaturated fatty acid, otherwise False.
        str: Reason for the classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject cyclic molecules (most fatty acids have an acyclic backbone)
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule is cyclic, which is not typical for a fatty acid chain"

    # Identify a carboxylic acid group.
    # This SMARTS pattern will match both protonated and deprotonated acid forms.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Molecule does not contain a carboxylic acid group, thus not a fatty acid"

    # For simplicity, take the first match. The first atom in the match is the acid carbon.
    acid_match = acid_matches[0]
    acid_carbon = acid_match[0]
    
    # Check that the carboxylic acid's carbon is terminal: it should have exactly one carbon neighbor.
    acid_atom = mol.GetAtomWithIdx(acid_carbon)
    carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Acid group is not terminal; the acid carbon has an unexpected bonding pattern"

    # Require that there is a sufficiently long alkyl chain.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 8:
        return False, "Molecule does not have enough carbon atoms to be considered a fatty acid"

    # Count carbon–carbon multiple bonds (double or triple) in the molecule,
    # ignoring the carbonyl C=O bond of the carboxyl and ensuring the unsaturation is internal.
    unsat_count = 0
    for bond in mol.GetBonds():
        # Only consider multiple bonds (double or triple)
        if bond.GetBondType() not in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        # Consider only bonds where both atoms are carbon.
        if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
            continue

        # Exclude the acid carbonyl bond (acid carbon to oxygen)
        # Note: Since the acid carbon is already identified and its carbonyl oxygen is not carbon,
        # this condition is mostly precautionary.
        if (a1.GetIdx() == acid_carbon or a2.GetIdx() == acid_carbon):
            # Skip if the other atom is not connected to further carbons (i.e. it's terminal)
            other_atom = a2 if a1.GetIdx() == acid_carbon else a1
            if other_atom.GetAtomicNum() != 6:
                continue

        # Ensure that the unsaturation is in the main chain:
        # Both atoms in the multiple bond should have at least one other carbon neighbor (i.e. not terminal methyl).
        a1_carbon_neighbors = [nbr for nbr in a1.GetNeighbors() if nbr.GetAtomicNum() == 6]
        a2_carbon_neighbors = [nbr for nbr in a2.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # If either end is terminal (only one carbon neighbor), skip counting this bond.
        if len(a1_carbon_neighbors) < 2 or len(a2_carbon_neighbors) < 2:
            continue

        unsat_count += 1

    if unsat_count != 1:
        return False, f"Fatty acid chain unsaturation count is {unsat_count} instead of the required 1"

    return True, ("Molecule contains a carboxylic acid group, is acyclic, has a sufficiently long fatty acid chain, "
                  "and exactly one internal C–C double/triple bond in the chain, classifying it as a monounsaturated fatty acid")