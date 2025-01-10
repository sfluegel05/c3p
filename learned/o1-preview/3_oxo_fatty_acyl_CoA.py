"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
"""
Classifies: 3-oxo-fatty acyl-CoA
Definition: An oxo fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any 3-oxo-fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Add hydrogens for accurate valence calculations
    mol_h = Chem.AddHs(mol)
    
    # Define the CoA substructure pattern (simplified)
    # This pattern includes the adenine moiety connected via ribose and phosphate groups
    # to the pantetheine moiety ending with the thiol group
    coa_smarts = """
    N1C=NC2=C1N=CN=C2N
    """
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not coa_pattern:
        return False, "Invalid CoA SMARTS pattern"

    # Check for CoA substructure
    if not mol_h.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # Find sulfur atoms connected to acyl carbonyl carbons (C(=O)-S linkage)
    for atom in mol_h.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Check for neighboring carbonyl carbon (C(=O))
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:  # Carbon atom
                    # Check if it's a carbonyl carbon (double-bonded to oxygen)
                    is_carbonyl = False
                    for bond in nbr.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            other_atom = bond.GetOtherAtom(nbr)
                            if other_atom.GetAtomicNum() == 8:  # Oxygen atom
                                is_carbonyl = True
                                break
                    if is_carbonyl:
                        acyl_carbon = nbr
                        # Now traverse to alpha carbon
                        alpha_carbons = [n for n in acyl_carbon.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
                        if len(alpha_carbons) != 1:
                            continue  # Not a typical fatty acyl chain
                        alpha_carbon = alpha_carbons[0]
                        # Traverse to beta carbon
                        beta_carbons = [n for n in alpha_carbon.GetNeighbors() if n.GetIdx() != acyl_carbon.GetIdx()]
                        if len(beta_carbons) != 1:
                            continue  # Not a typical fatty acyl chain
                        beta_carbon = beta_carbons[0]
                        # Check if beta carbon has a ketone group
                        has_ketone = False
                        for bond in beta_carbon.GetBonds():
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                other_atom = bond.GetOtherAtom(beta_carbon)
                                if other_atom.GetAtomicNum() == 8:  # Oxygen
                                    has_ketone = True
                                    break
                        if has_ketone:
                            return True, "Molecule is a 3-oxo-fatty acyl-CoA"
                        else:
                            return False, "Beta carbon does not have a ketone group"
    return False, "No matching 3-oxo-fatty acyl-CoA structure found"