"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol (A compound that contains exactly two free/alcoholic hydroxy groups)
Definition:
    “Diol” means that the molecule has exactly two free (alcoholic) –OH groups attached to sp3 carbons.
    –OH groups that are part of carbonyl-based functionalities (e.g. carboxylic acids) or attached to non-aliphatic carbons are ignored.
    
The algorithm:
  1. Parse the SMILES string and add explicit hydrogens.
  2. Use a SMARTS pattern "[OX2H]" to identify hydroxyl groups.
  3. For each –OH found, ensure that:
       a. It is attached to at least one carbon.
       b. The attached carbon is sp³ (aliphatic) and not double-bonded to another oxygen (i.e. not forming a carbonyl).
  4. Count the valid –OH groups; if exactly two are found, classify the molecule as a diol.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is classified as a diol (contains exactly two free/alcoholic hydroxyl groups)
    based on its SMILES string.
    
    The function:
      - Converts the SMILES string to an RDKit molecule and adds explicit hydrogens.
      - Uses the SMARTS pattern [OX2H] to find hydroxyl groups.
      - For each observed hydroxyl group, confirms that it is attached to an sp3 carbon,
        and that the carbon is not part of a carbonyl (i.e. not double-bonded to an oxygen).
      - Returns True if exactly two such groups are found.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a diol (exactly two qualifying –OH groups), otherwise False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens, ensuring the hydroxyl hydrogens are represented.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for –OH group (oxygen with 2 connections and an attached hydrogen)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    valid_oh_count = 0
    for match in oh_matches:
        oh_idx = match[0]
        oh_atom = mol.GetAtomWithIdx(oh_idx)
        # Get neighbors of the –OH oxygen
        neighbors = oh_atom.GetNeighbors()
        
        # The –OH must be attached to at least one carbon.
        # Additionally, we require that at least one attached carbon is sp3 (typical for aliphatic alcohols)
        has_aliphatic_carbon = False
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 6:
                # Check hybridization: only consider sp3 (non-aromatic/aliphatic)
                if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                    continue
                # Now check the bonds on the carbon to ensure it is not part of a carbonyl.
                # We skip the bond from this carbon to the –OH (our current oxygen).
                carbon_atom = nbr
                is_carbonyl = False
                for bond in carbon_atom.GetBonds():
                    # Get the other atom in the bond.
                    other = bond.GetOtherAtom(carbon_atom)
                    # Skip if the bond is the one connecting the hydroxyl oxygen.
                    if other.GetIdx() == oh_idx:
                        continue
                    # Look for a double bond from this carbon to an oxygen.
                    if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        is_carbonyl = True
                        break
                if not is_carbonyl:
                    has_aliphatic_carbon = True
                    break
        if has_aliphatic_carbon:
            valid_oh_count += 1
    
    if valid_oh_count == 2:
        return True, "Molecule contains exactly two free (alcoholic) hydroxyl groups attached to sp3 carbons and is classified as a diol."
    else:
        return False, f"Molecule contains {valid_oh_count} qualifying hydroxyl groups, which does not match the diol definition (exactly two required)."

# Example usage:
# test_smiles = "OCCCCCCCCCCCCO"  # Example: docosane-1,3-diol
# result, reason = is_diol(test_smiles)
# print(result, reason)