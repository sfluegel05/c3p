"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: The conjugate base of a fatty acid (i.e. fatty acid anion)
Definition: "The conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group of the corresponding fatty acid."
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    We require the molecule to possess a terminal deprotonated carboxylate group (i.e. [CX3](=O)[O-])
    and to have an acyl chain of sufficient length (we require at least three carbons in the chain besides the carboxylate carbon).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is identified as a fatty acid anion, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define the SMARTS for a deprotonated carboxylate group.
    # This pattern looks for a carbon (sp2) doubly bonded to an oxygen and singly bonded to an oxygen carrying negative charge.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No deprotonated carboxylate group found."

    # Get all matches of the carboxylate group.
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # We look for a candidate in which the carboxyl carbon (first atom of the match) is terminal,
    # i.e. it has only one carbon neighbor (its α–carbon). This is a heuristic to judge that it is coming
    # from a fatty acid rather than a carboxylate embedded in (for example) an amino acid.
    for match in matches:
        carboxylC = match[0]  # the carboxyl carbon as defined by SMARTS
        atom = mol.GetAtomWithIdx(carboxylC)
        # Find all neighbors that are carbon atoms (atomic number 6)
        carbon_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # We expect a terminal carboxylate: the carboxyl carbon should have exactly one carbon neighbor
        if len(carbon_neighbors) != 1:
            continue  # try next match if available

        # Now, verify that the alkyl (or acyl) chain is not too short.
        # We count all carbon atoms in the molecule but subtract the carboxyl carbon.
        total_carbons = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
        # The alkyl chain length is then roughly total carbons - 1.
        if (total_carbons - 1) < 3:
            # If there are less than 3 carbons in the chain, the fatty acyl chain is too short.
            return False, "Alkyl chain too short for a fatty acid."
        
        # Also, as a simple check for a fatty acid we might expect a predominance of carbon relative
        # to other heteroatoms. (This is only a heuristic.)
        heavy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1]  # ignore hydrogens
        carbon_ratio = total_carbons / len(heavy_atoms)
        if carbon_ratio < 0.4:
            return False, "Too many heteroatoms relative to carbons to be a fatty acid anion."

        # If these conditions are satisfied, we classify the molecule as a fatty acid anion.
        return True, "Molecule contains a terminal deprotonated carboxylate group with a sufficiently long acyl chain."

    # If none of the carboxylate matches satisfies the terminal condition, then it does not qualify.
    return False, "Carboxylate group found but not in a terminal (fatty acid-like) arrangement."

# Example usage (uncomment for testing):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCC([O-])=O"  # cerotate example
# print(is_fatty_acid_anion(test_smiles))