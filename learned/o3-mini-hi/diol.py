"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol (A compound that contains two hydroxy groups)
Definition (for this classifier):
    “Diol” means that the molecule has exactly two free (alcoholic) –OH groups.
    Here, we ignore –OH groups that are part of carboxylic acid (or similar carbonyl-based)
    functionalities. In addition, we require that the –OH groups are attached to a carbon.
    
Note: This classification is based on the idea that for a “diol” the two –OH groups should be
    alcoholic (i.e. not acidic like in a carboxyl group) and be attached to carbon atoms.
    The molecule is first “saturated” with explicit hydrogens to ensure that –OH groups are not missed.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is classified as a diol (contains exactly two free/alcoholic hydroxy groups)
    based on its SMILES string.
    
    The algorithm:
      1. Parses the SMILES and adds explicit hydrogens.
      2. Iterates over all oxygen atoms and selects those that:
           - Have at least one explicit hydrogen attached (indicating an –OH group)
           - Are attached to at least one carbon atom.
           - Are not attached (via that carbon) to a carbonyl group (i.e. the carbon is not double-bonded to another O)
             which would suggest a carboxylic acid functionality.
      3. If exactly two such “alcoholic” –OH groups are found, the molecule is classified as a diol.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a diol (exactly two qualifying –OH groups), False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES into an RDKit molecule and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    def is_alcoholic_oh(oh_atom):
        """
        Determines whether a given oxygen atom (assumed to be O) represents a free/alcoholic –OH group.
        It must have at least one hydrogen and be attached to at least one carbon.
        In addition, if the attached carbon shows a carbonyl double bond to another oxygen, then this –OH is likely
        part of a carboxylic acid and should not be counted.
        """
        # Check that the atom is oxygen
        if oh_atom.GetAtomicNum() != 8:
            return False
        # Check if it has at least one hydrogen attached.
        if oh_atom.GetTotalNumHs() < 1:
            return False
        
        # Look for a neighboring carbon atom.
        attached_carbons = [nbr for nbr in oh_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not attached_carbons:
            return False  # Not attached to any carbon
        
        # For each carbon neighbor, check if that carbon is involved in a C=O double bond.
        # If any attached carbon has a double-bonded oxygen (other than the current –OH), then skip this –OH.
        for carbon in attached_carbons:
            for bond in carbon.GetBonds():
                # Look for a bond (other than the one connecting to our –OH oxygen)
                other_atom = bond.GetOtherAtom(carbon)
                # Check for a double bond to another oxygen.
                if other_atom.GetAtomicNum() == 8 and other_atom.GetIdx() != oh_atom.GetIdx():
                    if bond.GetBondTypeAsDouble() == 2.0 or bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        # Found a carbonyl bond; treat this –OH as part of a carboxylic acid.
                        return False
        return True

    # Count the qualifying alcoholic hydroxy groups.
    alcoholic_oh_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            if is_alcoholic_oh(atom):
                alcoholic_oh_count += 1

    if alcoholic_oh_count == 2:
        return True, "Molecule contains exactly two free (alcoholic) hydroxy groups and is classified as a diol."
    else:
        return False, f"Molecule contains {alcoholic_oh_count} qualifying hydroxy groups, which does not match the diol definition (exactly two required)."

# Example usage (you can uncomment the following lines to test):
# test_smiles = "OC(C(O)(C)C)(C)C"  # pinacol
# result, reason = is_diol(test_smiles)
# print(result, reason)