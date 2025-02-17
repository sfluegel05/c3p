"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: Monounsaturated Fatty Acid

Definition: Any fatty acid with one double or triple bond in the fatty acid chain 
(and all remaining carbon–carbon bonds are single). Fatty acids are characterized 
by having a carboxylic acid group, and the monounsaturated ones have exactly one 
non–acid-group carbon–carbon unsaturation (double or triple bond).
"""

from rdkit import Chem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    
    The function performs the following steps:
      1. Parses the SMILES string.
      2. Checks for a carboxylic acid group (using a SMARTS pattern).
      3. Counts the unsaturation (double or triple C–C bonds) in the molecule, 
         ignoring the carbonyl double bond in the carboxylic acid group.
         
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

    # Identify a carboxylic acid group.
    # This SMARTS pattern matches a generic carboxylic acid; it matches both protonated and deprotonated forms.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Molecule does not contain a carboxylic acid group, thus not a fatty acid"

    # For simplicity, we will take the first carboxyl match.
    # In the match, the first atom is the carboxyl carbon.
    acid_carbon_idx = acid_matches[0][0]

    # Now count the carbon–carbon multiple bonds (double or triple) in the molecule.
    # We wish to ignore the C=O bond that belongs to the carboxyl group.
    unsat_count = 0
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type not in [Chem.BondType.DOUBLE, Chem.BondType.TRIPLE]:
            continue  # Only interested in multiple bonds

        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()

        # Consider only bonds where both atoms are carbon.
        if begin_atom.GetAtomicNum() != 6 or end_atom.GetAtomicNum() != 6:
            continue

        # If one of the atoms is the carboxyl carbon then check if the other atom is oxygen.
        # In the carboxyl group the carbonyl (C=O) bond is not considered as part of the fatty acid chain unsaturation.
        if (begin_atom.GetIdx() == acid_carbon_idx and end_atom.GetAtomicNum() == 8) or \
           (end_atom.GetIdx() == acid_carbon_idx and begin_atom.GetAtomicNum() == 8):
            continue

        # Count this bond as an unsaturation in the fatty acid chain.
        unsat_count += 1

    if unsat_count != 1:
        return False, f"Fatty acid chain unsaturation count is {unsat_count} instead of the required 1"

    return True, "Molecule contains a carboxylic acid group and exactly one C–C double/triple bond in the chain, classifying it as a monounsaturated fatty acid"